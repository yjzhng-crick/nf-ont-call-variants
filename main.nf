#!/usr/bin/env nextflow

/*
========================================================================================
   Variant calling Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-ont-call-variants
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println """\
         S C B I R   N A N O P O R E  V A R I A N T   C A L L I N G   P I P E L I N E
         ============================================================================
         Nextflow pipeline to call variants from Nanopore FASTQ files from bacterial clones 
         relative to a wildtype control.

         Usage:
            nextflow run scbirlab/nf-ont-call-variants --sample_sheet <csv> --inputs <dir> --genome_fasta <genome.fasta> --snpeff_database <organism-id>
            nextflow run scbirlab/nf-ont-call-variants -c <config-file>

         Required parameters:
            sample_sheet      Path to a CSV with information about the samples 
                                 and FASTQ files to be processed
            genome_fasta      Path to reference genome FASTA
            snpeff_database   Name of snpEff database to use for annotation.
                                 This should be derived from the same assembly as 
                                 `genome_fasta`. You can get a list of databases 
                                 using `java -jar snpEff database`. Database names 
                                 often end in the assembly name, such as `gca_000015005`, 
                                 which you can check matches your `genome_fasta`.

         Optional parameters (with defaults):  
            inputs             Directory containing inputs. Default: "./inputs". 
            gatk_image         Path to GATK executable.
            snpeff_url         URL to download snpEff.
            trim_qual = 10     For `cutadapt`, the minimum Phred score for trimming 3' calls
            min_length = 10    For `cutadapt`, the minimum trimmed length of a read. Shorter 
                                 reads will be discarded.

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}
if ( !params.genome_fasta ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to genome_fasta")
}
if ( !params.snpeff_database ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a snpeff_database")
}

working_dir = params.outputs

processed_o = "${working_dir}/processed"
tables_o = "${working_dir}/tables"
multiqc_o = "${working_dir}/multi_qc"

log.info """\
         S C B I R   V A R I A N T   C A L L I N G   P I P E L I N E
         ===========================================================
         software
            GATK           : ${params.gatk_image}
            Clair3         : ${params.clair3_image}
            snpEff         : ${params.snpeff_url}
         inputs
            input dir.     : ${params.inputs}
            sample sheet   : ${params.sample_sheet}
         genome locations
            FASTA          : ${params.genome_fasta}
            SnpEff organism: ${params.snpeff_database}
         trimming 
            quality        : ${params.trim_qual}
            minimum length : ${params.min_length}
         output            
            Processed      : ${processed_o}
            Tables         : ${tables_o}
            MultiQC        : ${multiqc_o}
         """
         .stripIndent()

dirs_to_make = [processed_o, tables_o, multiqc_o]

println  """
            Making directories: 
          """.stripIndent()

dirs_to_make.each { 
   println "${it}: " 
   println file(it).mkdirs() ? "OK" : "Cannot create directory: ${it}"
}


/*
========================================================================================
   Create Channels
========================================================================================
*/

csv_ch = Channel.fromPath( params.sample_sheet, 
                           checkIfExists: true )
                .splitCsv( header: true )

sample_ch = csv_ch.map { tuple( it.sample_id,
                                file( "${params.inputs}/${it.reads}", 
                                      checkIfExists: true ) ) }

genome_ch = Channel.fromPath( params.genome_fasta, 
                              checkIfExists: true )

clair3_ch = Channel.of( params.clair3_image )
gatk_ch = Channel.of( params.gatk_image )


/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   clair3_ch | LOAD_CLAIR3
   gatk_ch | LOAD_GATK

   sample_ch | FASTQC
   sample_ch | TRIM_CUTADAPT

   genome_ch | MINIMAP2_INDEX 
   TRIM_CUTADAPT.out.main.combine(MINIMAP2_INDEX.out) | MINIMAP2_ALIGN

   genome_ch | BOWTIE2_INDEX 
   TRIM_CUTADAPT.out.main.combine(BOWTIE2_INDEX.out) | BOWTIE2_ALIGN

   genome_ch | FAIDX 

   genome_ch.combine( FAIDX.out )
            .set { fasta_faidx }
   
   MINIMAP2_ALIGN.out.main | MOSDEPTH 
   MINIMAP2_ALIGN.out.main | SAMTOOLS_STATS

   MINIMAP2_ALIGN.out.main.combine( fasta_faidx )
                          .combine( LOAD_CLAIR3.out ) \
    | CLAIR3

   CLAIR3.out.main | NORMALIZE_VCFS

   NORMALIZE_VCFS.out.map { it[0] }.concat( NORMALIZE_VCFS.out.map { it[1] } ).collect() \
      | MERGE_VCFS 

   SNPEFF( MERGE_VCFS.out, 
           Channel.value( params.snpeff_url ), 
           Channel.value( params.snpeff_database ) )
   TO_TABLE( SNPEFF.out.main.combine( LOAD_GATK.out ),
             Channel.value( params.snpeff_url ) )

   TRIM_CUTADAPT.out.logs.concat(
         FASTQC.out.logs,
         MINIMAP2_ALIGN.out.logs,
         BOWTIE2_ALIGN.out.logs,
         SAMTOOLS_STATS.out.logs,
         SNPEFF.out.logs)
      .flatten()
      .unique()
      .collect() \
      | MULTIQC

}

process LOAD_GATK {

   label 'some_mem'

   input:
   val gatk_image

   output:
   path "*.sif"

   script:
   """
   singularity pull ${gatk_image}
   """
}

process LOAD_CLAIR3 {

   label 'some_mem'

   input:
   val clair3_image

   output:
   path "*.sif"

   script:
   """
   singularity pull ${clair3_image}
   """
}

/* 
 * Do quality control checks
 */
process FASTQC {

   label 'med_mem'

   tag "${sample_id}"

   input:
   tuple val( sample_id ), path ( reads )

   output:
   path( "*.zip" ), emit: logs

   script:
   """
   zcat ${reads} > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """
   stub:
   """
   zcat ${reads} | head -n1000 > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads ${task.cpus} ${sample_id}.fastq
   rm ${sample_id}.fastq
   """

}

/*
 * Trim adapters from reads
 */
process TRIM_CUTADAPT {

   tag "${sample_id}"

   publishDir( processed_o, 
               mode: 'copy',
               pattern: "*.log" )

   input:
   tuple val( sample_id ), path( reads )

   output:
   tuple val( sample_id ), path( "*.trimmed.fastq.gz" ), emit: main
   path "*.log" , emit: logs

   script:
   """
   ln -s ${reads} ${sample_id}.fastq.gz
   cutadapt \
		-q ${params.trim_qual},${params.trim_qual} \
		--report=full \
		-o ${sample_id}.trimmed.fastq.gz \
      --json=${sample_id}.cutadapt.json \
		${sample_id}.fastq.gz > ${sample_id}.cutadapt.log
   """
   stub:
   """
   zcat ${reads} | head -n1000 | gzip --best > ${sample_id}.fastq.gz
   cutadapt \
		-q ${params.trim_qual},${params.trim_qual} \
		--report=full \
		-o ${sample_id}.trimmed.fastq.gz \
      --json=${sample_id}.cutadapt.json \
		${sample_id}.fastq.gz > ${sample_id}.cutadapt.log
   """
}

/*
 * Index the reference genome for use by Bowtie2.
 */
process BOWTIE2_INDEX {
   
   tag "${genome}"
   
   input:
   path genome

   output:
   tuple val( "${genome.getBaseName()}" ), path( "*.bt2" )  

   script:
   """
   bowtie2-build ${genome} ${genome.getBaseName()}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process BOWTIE2_ALIGN {

   label 'med_mem'
   errorStrategy 'ignore'

   tag "${sample_id}"

   input:
   tuple val( sample_id ), path( reads ), val( genome_idx_base ), path( genome_idx )

   output:
   tuple val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   path "*.log", emit: logs

   script:
   """
   bowtie2 \
      -x ${genome_idx_base} \
      --rdg 10,10 --very-sensitive \
      --trim-to 3:2500 \
      -U ${reads} -S ${sample_id}.mapped.sam 2> ${sample_id}.bowtie2.log
   samtools sort ${sample_id}.mapped.sam \
      -O bam -l 9 -o ${sample_id}.mapped.bam
   samtools index ${sample_id}.mapped.bam
   rm ${sample_id}.mapped.sam
   """
}

/*
 * Get alignment stats
 */
process MOSDEPTH {

   tag "${sample_id}"

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile ), path( idx )

   output:
   path "*.txt", emit: logs
   path "*.html", emit: plots

   script:
   """
   mosdepth ${bamfile.getBaseName()} ${bamfile}
   git clone https://github.com/brentp/mosdepth.git && python mosdepth/scripts/plot-dist.py ${bamfile.getBaseName()}.mosdepth.global.dist.txt
   mv dist.html ${bamfile.getBaseName()}.mosdepth.global.dist.html
   rm -r mosdepth
   """
}

/*
 * Get alignment stats
 */
process SAMTOOLS_STATS {

   tag "${sample_id}"

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile ), path( idx )

   output:
   path "*.stats", emit: logs
   path "*.png", emit: plots

   script:
   """
   samtools stats ${bamfile} > ${bamfile.getBaseName()}.stats
   plot-bamstats -p ${bamfile.getBaseName()}_plot ${bamfile.getBaseName()}.stats
   """
}

/*
 * Index the reference genome for use by Minimap2.
 */
process MINIMAP2_INDEX {

   tag "${genome}"
   
   input:
   path genome

   output:
   path "*.mmi"

   script:
   """
   minimap2 -x map-ont -d ${genome.getBaseName()}-ont.mmi ${genome}
   """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process MINIMAP2_ALIGN {

   tag "${sample_id}" 

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( reads ), path( idx )

   output:
   tuple val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   path "*.log", emit: logs

   script:
   """
   minimap2 -aL --cs --M \
      -a ${idx} <(zcat ${reads}) \
      | samtools sort -@${task.cpus} -O bam -l 9 -o ${sample_id}.mapped.bam \
      2> ${sample_id}.minimap2.log
   samtools index ${sample_id}.mapped.bam
   """
}

/*
 * 
 */
process FAIDX {

   tag "${genome}"
   
   input:
   path genome

   output:
   path "*.fai" 

   script:
   """
   samtools faidx ${genome}
   """
}


process CLAIR3 {

   label 'big_cpu'

   tag "${sample_id}" 

   input:
   tuple val( sample_id ), path( bamfile ), path( idx ), path( fasta ), path( fai ), path( clair3_image )

   output:
   tuple val( sample_id ), path( "*.merge_output.vcf.gz" ), path( "*.merge_output.vcf.gz.tbi" ), path( fasta ), path( fai ), emit: main
   path "*.full_alignment.vcf.gz", emit: alignment
   path "*.log", emit: logs

   script:
   """
   THIS_DIR=\$(readlink -f \$(pwd))
   echo \$THIS_DIR
   git clone https://github.com/nanoporetech/rerio.git
   MODEL=r1041_e82_400bps_sup_v500
   python rerio/download_model.py --clair3 rerio/clair3_models/\$MODEL"_model"

   singularity exec \
      -B ${launchDir} \
      ${clair3_image} \
      /opt/bin/run_clair3.sh \
      --bam_fn=\$THIS_DIR/${bamfile} \
      --ref_fn=\$THIS_DIR/${fasta} \
      --sample_name=${sample_id}\
      --threads=${task.cpus} \
      --platform="ont" \
      --model_path="\$THIS_DIR/rerio/clair3_models/\$MODEL" \
      --output=\$THIS_DIR/${sample_id}_clair3 \
      --no_phasing_for_fa \
      --include_all_ctgs \
      --haploid_precise \
      --enable_long_indel

   mv ${sample_id}_clair3/merge_output.vcf.gz ${sample_id}_clair3.merge_output.vcf.gz
   mv ${sample_id}_clair3/merge_output.vcf.gz.tbi ${sample_id}_clair3.merge_output.vcf.gz.tbi
   mv ${sample_id}_clair3/full_alignment.vcf.gz ${sample_id}_clair3.full_alignment.vcf.gz
   mv ${sample_id}_clair3/run_clair3.log ${sample_id}_clair3.log
   rm -rf rerio
   """
}

process NORMALIZE_VCFS {

   label 'some_mem'

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( vcf ), path( vcf_idx ), path( fasta ), path( fai )

   output:
   tuple path( "*.normalized.vcf.gz" ), path( "*.tbi" )

   script:
   """
   zcat ${vcf} \
      | bcftools view -i 'GT="alt"' \
      | bcftools norm -f ${fasta} -a -c e -m - \
      | bcftools norm -aD \
      | bcftools +setGT - -- -t a -n c:M \
      | bcftools sort \
      | bcftools view -i 'GT="A"' -Oz -o ${sample_id}.normalized.vcf.gz
   bcftools index --tbi ${sample_id}.normalized.vcf.gz
   """

}


process MERGE_VCFS {

   label 'some_mem'

   input:
   path vcfs 

   output:
   path "merged.vcf"

   script:
   """
   bcftools merge --file-list <(ls -1 *.vcf.gz) -O v -o merged.vcf
   """
}


process SNPEFF {

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   path vcf
   val snpeff_url
   val snpeff_database

   output:
   path "*.ann.vcf", emit: main
   tuple path( "*.csv" ), path( "*.html" ), emit: logs

   script:
   """
   curl -v -L '${snpeff_url}' > snpEff_latest_core.zip
   unzip snpEff_latest_core.zip
   java -jar snpEff/snpEff.jar download -v ${snpeff_database}

   CHR_NAME=\$(grep -v '^#' ${vcf} | cut -f1 | sort -u)
   bcftools annotate \
      --rename-chrs <(printf \$CHR_NAME'\\tChromosome') \
      ${vcf} \
      > ${vcf}.1
   mv ${vcf}.1 ${vcf}

   java -jar snpEff/snpEff.jar -v -d \
      -ud 0 \
      -o gatk \
	   -csvStats ${vcf.getSimpleName()}.snpEff.csv \
      -stats  ${vcf.getSimpleName()}.snpEff.html \
      ${params.snpeff_database} ${vcf} \
      > ${vcf.getSimpleName()}.ann.vcf \
      && rm -r snpEff

   bcftools annotate \
      --rename-chrs <(printf 'Chromosome\\t'\$CHR_NAME) \
      ${vcf.getSimpleName()}.ann.vcf \
      > ${vcf.getSimpleName()}.ann.vcf.1
   mv ${vcf.getSimpleName()}.ann.vcf.1 ${vcf.getSimpleName()}.ann.vcf
   """
}

process TO_TABLE {

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple path( vcf ), path( gatk_image )
   val snpeff_url

   output:
   path "*.tsv"

   script:
   """
   curl -v -L '${snpeff_url}' > snpEff_latest_core.zip
   unzip snpEff_latest_core.zip
   java -jar snpEff/SnpSift.jar filter --inverse "(GEN[ALL].GT = GEN[0].GT)" ${vcf} \
      > ${vcf.getSimpleName()}.filtered.vcf

   java -jar snpEff/SnpSift.jar filter --inverse "(EFF[ALL].FUNCLASS = 'SILENT')" ${vcf.getSimpleName()}.filtered.vcf \
      > ${vcf.getSimpleName()}.high-eff.vcf

   set -x
   for f in ${vcf.getSimpleName()}.filtered.vcf ${vcf.getSimpleName()}.high-eff.vcf
   do
      OUTFILE=\$(basename \$f .vcf).tsv
      singularity exec \
         -B ${launchDir} \
         ${gatk_image} \
         /gatk/gatk VariantsToTable \
         -V \$f \
         -F CHROM -F POS -F REF -F ALT -F TYPE \
         -F QUAL -F FILTER -F EFF -GF GT \
         -O \$OUTFILE
      cat <(paste <(head -n1 \$OUTFILE) \
         <(printf 'Gene_Name\\tAmino_Acid_Change\\tFunctional Class')) \
         <(tail -n +2 \$OUTFILE | awk -F'|' 'BEGIN{OFS="\t"}{ print \$0,\$5,\$4,\$2 }') \
         > \$OUTFILE.1
      mv \$OUTFILE.1 \$OUTFILE
   done
   
   rm -r snpEff
   """
}

/*
 * Make log report
 */
process MULTIQC {

   publishDir( multiqc_o, 
               mode: 'copy' )

   input:
   path '*'

   output:
   tuple path( "*.html" ), path( "multiqc_data" )

   script:
   """
   multiqc .
   """
}


/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/