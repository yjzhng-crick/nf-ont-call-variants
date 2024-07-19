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
   genome_ch | PICARD_DICT 

   genome_ch.combine( FAIDX.out )
            .combine( PICARD_DICT.out )
            .set { fasta_faidx_dict }
   
   MINIMAP2_ALIGN.out.main | MOSDEPTH 
   MINIMAP2_ALIGN.out.main | SAMTOOLS_STATS

   MINIMAP2_ALIGN.out.main.combine( fasta_faidx_dict )
                          .combine( LOAD_CLAIR3.out ) \
    | CLAIR3

   CLAIR3.out.main.collect { it[1] }
                  .map { [it] }
                  .combine( CLAIR3.out.index.collect()
                                            .map { [it] } )
                  .combine( fasta_faidx_dict )
                  .combine( LOAD_GATK.out ) \
      | COMBINE_GENOTYPE_GVCFS | FILTER_VAR

   SNPEFF(FILTER_VAR.out, Channel.value(params.snpeff_url), Channel.value(params.snpeff_database))
   SNPEFF.out.main.combine( LOAD_GATK.out ) | TO_TABLE

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

   input:
   tuple val( sample_id ), path( reads ), path( idx )

   output:
   tuple val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   path "*.log", emit: logs

   script:
   """
   minimap2 -y --eqx -Y \
      --sam-hit-only --no-end-flt --sr --frag=yes \
      -F 3000 \
      -a ${idx} <(zcat ${reads}) \
      -o ${sample_id}.mapped.sam 2> ${sample_id}.minimap2.log
   samtools sort ${sample_id}.mapped.sam \
      -O bam -l 9 -o ${sample_id}.mapped.bam
   samtools index ${sample_id}.mapped.bam
   rm ${sample_id}.mapped.sam
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

/*
 * 
 */
process PICARD_DICT {

   tag "${genome}"
   
   input:
   path genome

   output:
   path "*.dict"

   script:
   """
   picard CreateSequenceDictionary \
      R=${genome} \
      O=${genome.getBaseName()}.dict
   """
}


process CLAIR3 {

   label 'big_cpu'

   tag "${sample_id}" 

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile ), path( idx ), path( fasta ), path( fai ), path( dict ), path( clair3_image )

   output:
   tuple val( sample_id ), path( "*.merge_output.gvcf.gz" ), path( "*.full_alignment.vcf.gz" ), emit: main
   path "*.merge_output.gvcf.gz.tbi", emit: index
   path "*.log", emit: logs

   script:
   """
   picard AddOrReplaceReadGroups \
      I=${bamfile} O=${sample_id}.rg0.bam \
      RGLB=${sample_id} RGPL="ONT" \
      RGPU=${sample_id} RGSM=${sample_id}
   samtools index ${sample_id}.rg0.bam

   THIS_DIR=\$(readlink -f \$(pwd))
   echo \$THIS_DIR
   singularity exec \
      -B ${launchDir} \
      ${clair3_image} \
      /opt/bin/run_clair3.sh \
      --bam_fn=\$THIS_DIR/${sample_id}.rg0.bam \
      --ref_fn=\$THIS_DIR/${fasta} \
      --sample_name=${sample_id}\
      --threads=${task.cpus} \
      --platform="ont" \
      --model_path="/opt/models/r941_prom_hac_g360+g422" \
      --output=\$THIS_DIR/${sample_id}_clair3 \
      --no_phasing_for_fa \
      --include_all_ctgs \
      --gvcf \
      --haploid_precise \
      && rm ${sample_id}.rg0.bam

   mv ${sample_id}_clair3/merge_output.gvcf.gz ${sample_id}_clair3.merge_output.gvcf.gz
   mv ${sample_id}_clair3/merge_output.gvcf.gz.tbi ${sample_id}_clair3.merge_output.gvcf.gz.tbi
   mv ${sample_id}_clair3/full_alignment.vcf.gz ${sample_id}_clair3.full_alignment.vcf.gz
   mv ${sample_id}_clair3/run_clair3.log ${sample_id}_clair3.log
   """
}

/*
 * GATK HTC
 */
process COMBINE_GENOTYPE_GVCFS {

   label 'some_mem'

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple path( gvcfs ), path( idxs ), path( fasta ), path( fai ), path( dict ), path( gatk_image )

   output:
   tuple path( "*.genotyped.vcf" ), path( fasta ), path( fai ), path( dict ), path( gatk_image )

   script:
   var gvcfs_ = gvcfs.each { it.getName() }.join(' --variant ')
   """
   singularity exec \
      -B ${launchDir} \
      ${gatk_image} \
      /gatk/gatk CombineGVCFs \
      -R ${fasta} --variant ${gvcfs_} -O combined.g.vcf
   singularity exec \
      -B ${launchDir} \
      ${gatk_image} \
      /gatk/gatk GenotypeGVCFs \
	   -R ${fasta} -V combined.g.vcf -O combined.genotyped.vcf
   """
}

/*
 * GATK HTC
 */
process FILTER_VAR {

   label 'some_mem'

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple path( vcf ), path( fasta ), path( fai ), path( dict ), path( gatk_image )

   output:
   path "*.vcf"

   script:
   """
   singularity exec \
      -B ${launchDir} \
      ${gatk_image} \
      /gatk/gatk VariantFiltration \
	   --filter-name "Auwera2013_QD" \
	   --filter-expression "QD<2.0" \
      --filter-name "Auwera2013_FS" \
		--filter-expression "FS>60.0" \
		--filter-name "Auwera2013_MQ" \
		--filter-expression "MQ<40.0" \
		--filter-name "Auwera2013_HaplotypeScore" \
		--filter-expression "HaplotypeScore>13.0" \
		--filter-name "Auwera2013_MappingQualityRankSum" \
		--filter-expression "MappingQualityRankSum<-12.5" \
		--filter-name "Auwera2013_ReadPosRankSum" \
		--filter-expression "ReadPosRankSum<-8.0" \
	   -R ${fasta} -V ${vcf} -O ${vcf.getSimpleName()}.filtered.vcf
   """
}

/*
 * GATK HTC
 */
process SNPEFF {

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   path gvcf
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

   CHR_NAME=\$(grep -v '^#' ${gvcf} | cut -f1 | sort -u)
   bcftools annotate \
      --rename-chrs <(printf \$CHR_NAME'\\tChromosome') \
      ${gvcf} \
      > ${gvcf}.1
   mv ${gvcf}.1 ${gvcf}

   java -jar snpEff/snpEff.jar -v -d \
      -o gatk \
      -ud 0 \
	   -csvStats ${gvcf.getSimpleName()}.snpEff.csv \
      -stats  ${gvcf.getSimpleName()}.snpEff.html \
      ${params.snpeff_database} ${gvcf} \
      > ${gvcf.getSimpleName()}.ann.vcf \
      && rm -r snpEff

   bcftools annotate \
      --rename-chrs <(printf 'Chromosome\\t'\$CHR_NAME) \
      ${gvcf.getSimpleName()}.ann.vcf \
      > ${gvcf.getSimpleName()}.ann.vcf.1
   mv ${gvcf.getSimpleName()}.ann.vcf.1 ${gvcf.getSimpleName()}.ann.vcf
   """
}

/*
 * GATK HTC
 */
process TO_TABLE {

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple path( vcf ), path( gatk_image )

   output:
   path "*.tsv"

   script:
   """
   singularity exec \
      -B ${launchDir} \
      ${gatk_image} \
      /gatk/gatk VariantsToTable \
     	-V ${vcf} \
     	-F CHROM -F POS -F REF -F ALT -F TYPE \
		-F QUAL -F FILTER -F EFF -GF GT \
     	-O ${vcf.getSimpleName()}0.tsv

   cat <(paste <(head -n1 ${vcf.getSimpleName()}0.tsv) <(printf 'Gene_Name\\tAmino_Acid_Change\\tFunctional Class')) \
       <(tail -n +2 ${vcf.getSimpleName()}0.tsv | awk -F'|' 'BEGIN{OFS="\t"}{ print \$0,\$5,\$4,\$2 }') \
       > ${vcf.getSimpleName()}.tsv

   #rm ${vcf.getSimpleName()}0.tsv

   NCOL=\$(( \$(head -n1 ${vcf.getSimpleName()}.tsv | awk -F'\\t' '{ print NF }') - 11 ))
   AWK_COMMAND='(\$(NF - 4) != \$(NF - 3)) { print }'
   
   for i in \$(seq 2 \$NCOL)
   do
      AWK_COMMAND='(\$(NF - '\$(( \$i + 3 ))') != \$(NF - 3)) && '\$AWK_COMMAND
   done

   echo "\$AWK_COMMAND"

   cat <(head -n1 ${vcf.getSimpleName()}.tsv) \
       <(tail -n +2 ${vcf.getSimpleName()}.tsv \
         | awk -F'\\t' "\$AWK_COMMAND") \
       >  ${vcf.getSimpleName()}.filtered.tsv

   cat <(head -n1 ${vcf.getSimpleName()}.filtered.tsv) \
       <(tail -n +2 ${vcf.getSimpleName()}.filtered.tsv \
         | grep -v '|SILENT|' | sort -k6 -n -r ) \
      > ${vcf.getSimpleName()}.high-eff.tsv
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