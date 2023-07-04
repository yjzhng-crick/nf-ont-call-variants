# Nanopore variant calling pipeline

Nextflow pipeline to call variants from Nanopore FASTQ files from bacterial clones relative to a wildtype control.

The pipeline broadly recapitualtes, where possible, the GATK best practices for germline short variant calling. 

## Processing steps

For each sample:

1. Quality Trim reads using `cutadapt`. 
2. Map to genome FASTA using `minimap2`.
3. Mark duplicates with `picard MarkDuplicates`.
4. Re-align in "active" regions and calculate variant likelihood with GATK `HaplotypeCaller`.

Then merge resulting GVCFs using GATK `CombineGVCFs`. With the combined variant calls:

1. Joint genoptype with GATK `GenotypeGVCFs`.
2. Filter variants using GATK `VariantFiltration`.
3. Annotate variant effects using `snpEff`.
4. Filter out variants where all samples are identical to the wildtype control, which [is assumed to be the `sample_id` which is alphabetically last](#sample-sheet).
5. Write to output TSV.

### Other steps

1. Get FASTQ quality metrics with `fastqc`.
2. Generate alignment statistics with `samtools stats`.
2. Map to genome FASTA using `bowtie2` because `minimap2` logs are not compatible with `multiqc`. This way, some kind of alignment metrics are possible.
3. Compile the logs of processing steps into an HTML report with `multiqc`.

## Requirements

### Software

- Nextflow
    - At Crick, activate using `module load Nextflow`
    - Otherwise, [see below](#first-time-using-nextflow)
- `conda` or `mamba`
    - If possible, use `mamba` because it will be faster.
    - At Crick, activate using `module load Anaconda3`
- GATK
    - Download [here](https://github.com/broadinstitute/gatk/releases), and provide the path as `--gatk_path` or in the `nextflow.config` file ([see below](#inputs))
- Picard
    - Download [here](https://broadinstitute.github.io/picard/), and provide the path as `--picard_path` or in the `nextflow.config` file ([see below](#inputs))
- snpEff
    - Download [here](https://pcingola.github.io/SnpEff/), and provide the path as `--snpeff_path` or in the `nextflow.config` file ([see below](#inputs))

### Reference genome

You also need the genome FASTA and GFF annotations for the bacteria you are sequencing. These can be obtained from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/):

1. Search for your strain of interest, and open its main page
2. On the right-hand side, click `Customize view`, then `Customize` and check `Show sequence`. Finally, click `Update view`. You may have to wait a few minute while the sequence downloads.
3. Click `Send to: > Complete record > File > FASTA > Create file`
4. Save the files to a path which you provide as `--genome_fasta` [below](#inputs).

### First time using Nextflow?

If it's your first time using Nextflow on your system, you can install it using `conda`:

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

## Quick start

Make sure you have [GATK, Picard, and snpEff on your system](#software), and provide their paths as [parameters on the command line or in your `nextflow.config` file](#inputs).

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
nextflow run scbirlab/nf-ont-call-variants -r v0.0.1
```

For help, use `nextflow run scbirlab/nf-ont-call-variants --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV with information about the samples and FASTQ files to be processed
- `gatk_path`: path to GATK executable
- `picard_path`: path to Picard executable
- `snpeff_path`: path to snpEff executable
- `genome_fasta`: path to reference genome FASTA
- `snpeff_database`: name of snpEff database to use for annotation. This should be derived from the same assembly as `genome_fasta`. You can get a list of databases using `java -jar snpEff database`. Database names often end in the assembly name, such as `gca_000015005`, which you can check matches your `genome_fasta`

The following parameters have default values which can be overridden if necessary.

- `trim_qual = 10` : For `cutadapt`, the minimum Phred score for trimming 3' calls
- `min_length = 10` : For `cutadapt`, the minimum trimmed length of a read. Shorter reads will be discarded

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
   
    gatk_path = "/path/to/gatk"
    picard_path = "/path/to/picard.jar"
    snpeff_path = "/path/to/snpEff.jar"

    sample_sheet = "/path/to/sample-sheet.csv"
    
    genome_fasta = "/path/to/MsmMC2155-CP000480.1.fasta"
    snpeff_database = "Mycolicibacterium_smegmatis_mc2_155_gca_000015005"

}
```

Alternatively, you can provide the parameters on the command line:

```bash
nextflow run scbirlab/nf-ont-call-variants \
    --sample_sheet /path/to/sample-sheet.csv \
    --gatk_path /path/to/gatk \
    --picard_path /path/to/picard.jar \
    --snpeff_path /path/to/snpEff.jar \
    --genome_fasta /path/to/MsmMC2155-CP000480.1.fasta \
    --snpeff_database Mycolicibacterium_smegmatis_mc2_155_gca_000015005
``` 

### Sample sheet

The sample sheet is a CSV file providing information about which FASTQ files belong to which sample.

The file must have a header with the column names below, and one line per sample to be processed.

- `sample_id`: the unique name of the sample. The wildtype must be **named so that it is alphabetically last**
- `reads`: path to compressed FASTQ files derived from Nanopore sequencing

Here is an example of the sample sheet:

| sample_id | reads                                  |
| --------- | -------------------------------------- | 
| wt        | /path/to/reads/WT/raw_reads.fastq.gz   | 
| mut1      | /path/to/reads/mut1/raw_reads.fastq.gz | 

## Outputs

Outputs are saved in the same directory as `sample_sheet`. They are organised under three directories:

- `processed`: FASTQ files and logs resulting from alignments
- `tables`: tables and VCF files corresponding to variant calls
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