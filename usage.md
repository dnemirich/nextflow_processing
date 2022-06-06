The pipeline includes several steps:

1. Preparation step.
2. Transcriptomes processing and genotyping.
3. Exomes processing, genotyping and annotation.

I. PREPARATION STEP

In order to start, launch **preparation.sh** script. This script will install nextflow and other tools, required for subsequent analysis. It will create conda environment called **melanoma** and make all necessary directories.


The recommended structure for input directory is as follows:
```sh
./input
└── sample1
    ├── exomes
    │   ├── normal
    │   │   ├── normal.1.fastq.gz
    │   │   └── normal.2.fastq.gz
    │   └── tumor
    │       ├── tumor.1.fastq.gz
    │       └── tumor.2.fastq.gz
    └── transcriptomes
        ├── filtered_transcriptomes
        │   ├── tr_tumor.1.1_filtered.noadapt.fastq.gz
        │   ├── tr_tumor.1.2_filtered.noadapt.fastq.gz
        │   ├── tr_tumor.2.1_filtered.noadapt.fastq.gz
        │   └── tr_tumor.2.2_filtered.noadapt.fastq.gz
        ├── tr_tumor.1.1.fastq.gz
        ├── tr_tumor.1.2.fastq.gz
        ├── tr_tumor.2.1.fastq.gz
        └── tr_tumor.2.2.fastq.gz
```

The reference directory should contain **hg38.analysisSet.fa** file. The data directory should contain **gencode.v31.chr_patch_hapl_scaff.annotation.gtf** and **ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf** files. To process aforementioned files:
```sh
conda activate melanoma
bash wgs_preparation.sh 
```

Transcriptomes filtering is a separate step. In order to perform filtering, launch **reads_filtering.sh** script. The following arguments are required:

```sh
conda activate melanoma #in case you have not activated it
bash reads_filtering.sh -f <r1.fastq.gz> -r <r2.fastq.gz> -T <FastqFilterer threshold value> -t <tiles to remove, comma separated> \
-a <first adapter to remove> -A <second adapter to remove> -c <number of cores> -o <path to output directory>
```

II. TRANSCRIPTOMES PROCESSING.

Transcriptomes can be processed with **human_nf_PART1.nf** script. You can specify parameters either through direct edit of the script (INPUT FILES and OUTPUT CONFIGURATION sections) or through usage of appropriate flags.

In case you edit the script directly, you should run pipeline this way:

```sh

nextflow run human_nf_PART1.nf

```


Otherwise, you need to specify arguments:

```sh

nextflow run human_nf_PART1.nf --filtered '<path_to_filtered_transcriptomes/*{1,2}_filtered.noadapt.fastq.gz>' --qc <output path to qc reports>\
--sortDir <output path to directory for sorted aligned and unaligned transcriptomes> \
--sm <sample_name in the following format **sample1_tumor_rna**> --merged <output path to softly-clipped transcriptomes> \
--genotyping <output path to results of genotyping> --outdir <path to directory for all intermediate files>

```

**!!!! You do not need to create output directories by yourself.**

Example:

```sh

nextflow run human_nf_PART1.nf --filtered 'input/sample1/transcriptomes/filtered_transcriptomes/*{1,2}_filtered.noadapt.fastq.gz' \
 --qc qc/sample1 --sortDir sorted/sample1 --sm sample1_tumor_rna --merged merged/sample1 --genotyping genotyping/sample1 --outdir intermediate
```

III. EXOMES PROCESSING.


Exomes processing and annotation can be performed with **human_nf_PART2.nf** script. The idea is the same: you can either edit the input directly or provide necessary arguments with flags.

```sh
nextflow run human_nf_PART2.nf --normal '<path_to_exomes_for_normal_tissue/*{1,2}.fastq.gz>' \
--tumor '<path_to_exomes_for_cancer_tissue>/*{1,2}.fastq.gz' --sm <sample name in the following format **sample1_**> \
--segments <path to tsv file with tumor segments> --qc <output path to qc reports directory> \
--outdir <path to directory for all internediate files> --aligned <output path to directory with aligned exomes> \
--geno <output path to results of genotyping>
```
Example:

```sh
nextflow run human_nf_PART2.nf --normal 'input/sample1/exomes/normal/*{1,2}.fastq.gz' \
--tumor 'input/sample1/exomes/tumor/*{1,2}.fastq.gz' --sm sample1_ \
--segments genotyping/sample1/segments.tsv --qc qc/sample1 --outdir intermediate_files/sample1 \
--aligned aligned/sample1 --geno genotyping/sample1

```

