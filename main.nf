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
         S C B I R   N A N O P O R E  s b R N A - S E Q   P I P E L I N E
         ================================================================
         Nextflow pipeline to process raw Nanopore POD5 files from 
         multiple bacterial samples into a gene x cell count table. Applies 
         duplex basecalling, therefore requires v10.4 chemistry.

         Usage:
            nextflow run sbcirlab/nf-ont-sbrnaseq --sample_sheet <csv> --data_dir <dir> --genome_fasta_dir <dir> --genome_gff_dir <dir> --guppy_path <path>
            nextflow run sbcirlab/nf-ont-sbrnaseq -c <config-file>

         Required parameters:
            sample_sheet      Path to a CSV with information about the samples 
                                 and FASTQ files to be processed
            gatk_path         Path to GATK executable.
            picard_path       Path to Picard executable
            snpeff_path       Path to snpEff executable
            genome_fasta      Path to reference genome FASTA
            snpeff_database   Name of snpEff database to use for annotation.
                                 This should be derived from the same assembly as 
                                 `genome_fasta`. You can get a list of databases 
                                 using `java -jar snpEff database`. Database names 
                                 often end in the assembly name, such as `gca_000015005`, 
                                 which you can check matches your `genome_fasta`.

         Optional parameters (with defaults):   
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
if ( !params.gatk_path ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a gatk_path")
}
if ( !params.picard_path ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a picard_path")
}
if ( !params.snpeff_path ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a snpeff_path")
}
if ( !params.genome_fasta ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to genome_fasta")
}
if ( !params.snpeff_database ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a snpeff_database")
}

wd = file(params.sample_sheet)
working_dir = wd.getParent()

processed_o = "${working_dir}/processed"
tables_o = "${working_dir}/tables"
multiqc_o = "${working_dir}/multi_qc"

log.info """\
         S C B I R   V A R I A N T   C A L L I N G   P I P E L I N E
         ===========================================================
         software
            GATK           : ${params.gatk_path}
            Picard         : ${params.picard_path}
            snpEff         : ${params.snpeff_path}
         inputs
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
                                file( it.reads ) ) }

genome_ch = Channel.fromPath( params.genome_fasta, 
                              checkIfExists: true )

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   sample_ch | FASTQC
   sample_ch | TRIM_CUTADAPT

   genome_ch | MINIMAP2_INDEX 
   TRIM_CUTADAPT.out.main.combine(MINIMAP2_INDEX.out) | MINIMAP2_ALIGN

   genome_ch | BOWTIE2_INDEX 
   TRIM_CUTADAPT.out.main.combine(BOWTIE2_INDEX.out) | BOWTIE2_ALIGN

   genome_ch | FAIDX 
   genome_ch | PICARD_DICT 

   genome_ch.combine(FAIDX.out)\
            .combine(PICARD_DICT.out) \
            .set { fasta_faidx_dict }
   
   MINIMAP2_ALIGN.out.main | MARK_DUPES 
   MARK_DUPES.out.main | SAMTOOLS_STATS

   MARK_DUPES.out.main.combine(fasta_faidx_dict) \
      | HAPLOTYPE_CALLER

   HAPLOTYPE_CALLER.out.collect { it[1] } \
                   .map { [it] }
                   .combine(fasta_faidx_dict) \
      | COMBINE_GENOTYPE_GVCFS | FILTER_VAR

   FILTER_VAR.out | SNPEFF
   SNPEFF.out.main | TO_TABLE

   TRIM_CUTADAPT.out.logs.concat(
         FASTQC.out.logs,
         MINIMAP2_ALIGN.out.logs,
         BOWTIE2_ALIGN.out.logs,
         MARK_DUPES.out.logs,
         SAMTOOLS_STATS.out.logs,
         SNPEFF.out.logs)
      .flatten()
      .unique()
      .collect() \
      | MULTIQC

}

/* 
 * Do quality control checks
 */
process FASTQC {

   memory '32GB'

   tag "${sample_id}"

   input:
   tuple val( sample_id ), path ( reads )

   output:
   path( "*.zip" ), emit: logs

   script:
   """
   zcat ${reads} > ${sample_id}.fastq
   fastqc --noextract --memory 10000 --threads 4 ${sample_id}.fastq
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
               pattern: "*.{log,json}" )

   input:
   tuple val( sample_id ), path( reads )

   output:
   tuple val( sample_id ), path( "*.trimmed.fastq.gz" ), emit: main
   tuple path( "*.log" ),  path( "*.json" ), emit: logs

   script:
   """
   ln -s ${reads} ${sample_id}.fastq.qz
   cutadapt \
		-q ${params.trim_qual},${params.trim_qual} \
		--report=full \
		-o ${sample_id}.trimmed.fastq.gz \
      --json=${sample_id}.cutadapt.json \
		${sample_id}.fastq.qz > ${sample_id}.cutadapt.log
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

   memory '32GB'
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
process SAMTOOLS_STATS {

   tag "${sample_id}"

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile ), val( idx )

   output:
   path "*.stats", emit: logs

   script:
   """
   samtools stats ${bamfile} > ${bamfile.getBaseName()}.stats
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
   path( "*.mmi" )

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
 * Picard Mark Duplicates
 */
process MARK_DUPES {

   tag "${sample_id}" 

   publishDir( path: processed_o, 
               mode: 'copy', 
               pattern: "*.marked-dupes.bam" )

   input:
   tuple val( sample_id ), path( bamfile ), path( idx )

   output:
   tuple val( sample_id ), path( "*.bam" ), path( "*.bai" ), emit: main
   path "*.txt", emit: logs

   script:
   """
   picard MarkDuplicates \
      I=${bamfile} \
      O=${bamfile.getSimpleName()}.marked-dupes.bam \
      M=${bamfile.getSimpleName()}.marked-dupes.txt
   samtools index ${bamfile.getSimpleName()}.marked-dupes.bam
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
   path( "*.fai" )

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
   path( "*.dict" )

   script:
   """
   picard CreateSequenceDictionary \
      R=${genome} \
      O=${genome.getBaseName()}.dict
   """
}

/*
 * GATK HTC
 */
process HAPLOTYPE_CALLER {

   tag "${sample_id}" 

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple val( sample_id ), path( bamfile ), path( idx ), path( fasta ), path( fai ), path( dict )

   output:
   tuple val( sample_id ), path( "*.g.vcf" ), path( "*.htc.bam" )

   script:
   """
   picard AddOrReplaceReadGroups \
      I=${bamfile} O=${sample_id}.rg0.bam \
      RGLB=${sample_id} RGPL="ONT" \
      RGPU=${sample_id} RGSM=${sample_id}

   # hack to deal with tp:A:P tag from minimap2 not compatible with GATK
   samtools view -h ${sample_id}.rg0.bam | cut -f23 --complement | samtools sort - -l9 -o ${sample_id}.rg.bam
   samtools index ${sample_id}.rg.bam

   ${params.gatk_path} HaplotypeCaller \
      -R ${fasta} \
      -ERC GVCF \
      -ploidy 1 \
      --pcr-indel-model NONE \
      -A AlleleFraction \
      -A BaseQuality \
      -I ${sample_id}.rg.bam \
      -O ${sample_id}.g.vcf \
      -bamout ${sample_id}.htc.bam

   rm ${sample_id}.rg0.bam
   """
}

/*
 * GATK HTC
 */
process COMBINE_GENOTYPE_GVCFS {

   memory '16GB'

   publishDir( path: processed_o, 
               mode: 'copy' )

   input:
   tuple path( gvcfs ), path( fasta ), path( fai ), path( dict )

   output:
   tuple path( "*.genotyped.vcf" ), path( fasta ), path( fai ), path( dict )

   script:
   var gvcfs_ = gvcfs.each { it.getName() }.join(' --variant ')
   """
   ${params.gatk_path} CombineGVCFs \
	    -R ${fasta} --variant ${gvcfs_} -O combined.g.vcf
   ${params.gatk_path} GenotypeGVCFs \
	    -R ${fasta} -V combined.g.vcf -O combined.genotyped.vcf
   """
}

/*
 * GATK HTC
 */
process FILTER_VAR {

   memory '16GB'

   publishDir( path: tables_o, 
               mode: 'copy' )

   input:
   tuple path( vcf ), path( fasta ), path( fai ), path( dict )

   output:
   path( "*.vcf" )

   script:
   """
   ${params.gatk_path} \
      VariantFiltration \
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
   path( gvcf )

   output:
   path( "*.ann.vcf" ), emit: main
   tuple path( "*.csv" ), path( "*.html" ), emit: logs

   script:
   """
   CHR_NAME=\$(grep -v '^#' ${gvcf} | cut -f1 | sort -u)
   bcftools annotate \
      --rename-chrs <(printf \$CHR_NAME'\\tChromosome') \
      ${gvcf} \
      > ${gvcf}.1
   mv ${gvcf}.1 ${gvcf}

   java -jar ${params.snpeff_path} -v -d \
      -o gatk \
      -ud 0 \
	   -csvStats ${gvcf.getSimpleName()}.snpEff.csv \
      -stats  ${gvcf.getSimpleName()}.snpEff.html \
      ${params.snpeff_database} ${gvcf} \
      > ${gvcf.getSimpleName()}.ann.vcf

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
   path( vcf )

   output:
   path( "*.tsv" )

   script:
   """
   ${params.gatk_path} VariantsToTable \
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