# Nanopore variant calling pipeline

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-ont-call-variants/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

**scbirlab/nf-ont-call-variants** is a Nextflow pipeline to call variants from Nanopore FASTQ files from bacterial clones relative to a wildtype control.

The pipeline broadly recapitualtes, where possible, the GATK best practices for germline short variant calling. 

**Table of contents**

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

## Processing steps

For each sample:

1. Quality Trim reads using `cutadapt`. 
2. Map to genome FASTA using `minimap2`.
3. Mark duplicates with `picard MarkDuplicates`.
4. Call variants with `Clair3`.

Then merge resulting GVCFs using GATK `CombineGVCFs`. With the combined variant calls:

1. Joint genoptype with GATK `GenotypeGVCFs`.
2. Filter variants using GATK `VariantFiltration`.
3. Annotate variant effects using `snpEff`.
4. Filter out variants where all samples are identical to the wildtype control, which [is assumed to be the `sample_id` which is alphabetically last](#sample-sheet).
5. Write to output TSV.

### Other steps

1. Get FASTQ quality metrics with `fastqc`.
2. Generate alignment statistics and plots with `samtools stats`.
2. Map to genome FASTA using `bowtie2` because `minimap2` logs are not compatible with `multiqc`. This way, some kind of alignment metrics are possible.
3. Compile the logs of processing steps into an HTML report with `multiqc`.

## Requirements

### Software

You need to have Nextflow and `conda` installed on your system.

#### First time using Nextflow?

If you're at the Crick or your shared cluster has it already installed, try:

```bash
module load Nextflow
```

Otherwise, if it's your first time using Nextflow on your system, you can install it using `conda`:

```bash
conda install -c bioconda nextflow 
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

### Reference genome

You also need the genome FASTA for the bacteria you are sequencing. These can be obtained from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/):

1. Search for your strain of interest, and open its main page
2. On the right-hand side, click `Customize view`, then `Customize` and check `Show sequence`. Finally, click `Update view`. You may have to wait a few minute while the sequence downloads.
3. Click `Send to: > Complete record > File > FASTA > Create file`
4. Save the files to a path which you provide as `--genome_fasta` [below](#inputs).

## Quick start

Make a [sample sheet (see below)](#sample-sheet) and, optionally, a [`nextflow.config` file](#inputs) in the directory where you want the pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-ont-call-variants
```

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version which will not be automatically updated. If you want to ensure that you're using the very latest version of the pipeline, use the `-latest` flag.

```bash 
nextflow run scbirlab/nf-ont-call-variants -latest
```
If you want to run a particular tagged version of the pipeline, such as `v0.0.1`, you can do so using

```bash 
nextflow run scbirlab/nf-ont-call-variants -r v0.0.2
```

For help, use `nextflow run scbirlab/nf-ont-call-variants --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV with information about the samples and FASTQ files to be processed
- `genome_fasta`: path to reference genome FASTA
- `snpeff_database`: name of snpEff database to use for annotation. This should be derived from the same assembly as `genome_fasta`. Database names often end in the assembly name, such as `gca_000015005`, which you can check matches that of your `genome_fasta`

The following parameters have default values which can be overridden if necessary.

- `inputs = "inputs"` : The folder containing your inputs.
- `trim_qual = 10` : For `cutadapt`, the minimum Phred score for trimming 3' calls
- `min_length = 10` : For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded
- `gatk_image = "docker://broadinstitute/gatk:latest"` : Which GATK4 image to use
- `snpeff_url = "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"` : Where to download snpEff from
- `clair3_image = "docker://hkubal/clair3:latest"` : Which Clair3 image to use

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
    sample_sheet = "/path/to/sample-sheet.csv"
    inputs = "/path/to/inputs"
    genome_fasta = "/path/to/1_ASM28329v1_genomic.fna"
    snpeff_database = "Mycobacterium_smegmatis_str_mc2_155_gca_000283295"

}
```

Alternatively, you can provide the parameters on the command line:

```bash
nextflow run scbirlab/nf-ont-call-variants \
    --sample_sheet /path/to/sample-sheet.csv \
    --inputs /path/to/inputs
    --genome_fasta /path/to/1_ASM28329v1_genomic.fna \
    --snpeff_database Mycobacterium_smegmatis_str_mc2_155_gca_000283295
``` 

### Sample sheet

The sample sheet is a CSV file providing information about which FASTQ files belong to which sample.

The file must have a header with the column names below, and one line per sample to be processed.

- `sample_id`: the unique name of the sample. The wildtype must be **named so that it is alphabetically last**
- `reads`: path (relative to `inputs` option above) to compressed FASTQ files derived from Nanopore sequencing

Here is an example of the sample sheet:

| sample_id | reads                |
| --------- | -------------------- | 
| wt        | raw_reads.fastq.gz   | 
| mut1      | raw_reads.fastq.gz   | 

## Outputs

Outputs are saved in the same directory as `sample_sheet`. They are organised under three directories:

- `processed`: FASTQ files and logs resulting from alignments
- `tables`: tables, plots, and VCF files corresponding to variant calls
- `multiqc`: HTML report on processing steps

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-ont-call-variants/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [multiqc](https://multiqc.info/)
- [nextflow](https://www.nextflow.io/docs/latest/index.html)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
- [minimap2](https://lh3.github.io/minimap2/minimap2.html)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools](http://www.htslib.org/doc/samtools.html)
- [GATK](https://gatk.broadinstitute.org/hc/en-us)
- [Picard](https://broadinstitute.github.io/picard/)
- [snpEff](https://pcingola.github.io/SnpEff/)