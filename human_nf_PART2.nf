
// INPUT FILES

params.normal = "$baseDir/input/sample1/exomes/normal/*{1,2}.fastq.gz"
params.tumor = "$baseDir/input/sample1/exomes/tumor/*{1,2}.fastq.gz"
params.reference = "$baseDir/reference/hg38.analysisSet.fa"
params.referenceDir = "$baseDir/reference"
params.dataDir = "$baseDir/data"
params.sm = "sample1_"
params.segments = "$baseDir/genotyping/sample1/segments.tsv"
params.wgs = "$baseDir/data/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf"
params.funcotator_data = "$baseDir/data/funcotator/funcotator_*"


// OUTPUT CONFIGURATION

params.qc = "$baseDir/qc/sample1"
params.outdir = './intermediate_files/sample1'
params.aligned = "$baseDir/aligned/sample1"
params.geno = "$baseDir/genotyping/sample1"


Channel
 .fromFilePairs( params.normal )
 .ifEmpty { error "Cannot find any reads matching: ${params.normal}" }
 .into { normal_pairs_ch; normal_pairs_ch_2}

Channel
 .fromFilePairs(params.normal, flat: true)
 .set { separate_normal }

 Channel
  .fromFilePairs( params.tumor )
  .ifEmpty { error "Cannot find any reads matching: ${params.tumor}" }
  .into { tumor_pairs_ch; tumor_pairs_ch_2}

  Channel
   .fromFilePairs(params.tumor, flat: true)
   .set { separate_tumor }



 /*
  * QC STEP
  */

  process fastqc{

    publishDir "${params.qc}", mode: 'copy'

    input:
    set val(id), file(reads) from normal_pairs_ch
    set val(id_2), file(reads_2) from tumor_pairs_ch


    output:
    file "*.{zip,html,txt}" into fastqc_results

    script:

    """
    fastqc -t 10 $reads $reads_2
    """
  }


  /*
   * REFERENCE PREPARATION STEP
   */


   process reference_preparation_1{

     tag "$params.reference"
     publishDir "${params.refDir}", mode: 'copy'

     input:
     path reference from params.reference

     output:
     file "*.fai" into reference_index_ch

     script:
     """
     samtools faidx $reference
     """
   }


   process reference_preparation_2{

     tag "$params.reference"
     publishDir "${params.refDir}", mode: 'copy'

     input:
     path reference from params.reference

     output:
     path "*.dict" into reference_dict_ch

     script:
     """
     picard CreateSequenceDictionary R= $reference O= ${reference.baseName}.dict
     """
   }


  process reference_preparation_3{

    tag "$params.reference"
    publishDir "${params.refDir}", mode: 'copy'

    input:
    path reference from params.reference

    output:
    file "*.{amb,ann,bwt,pac,sa}" into reference_index_ch_2

    script:
    """
    bwa index $reference
    """
  }


  /*
  *  EXOMES PROCESSING
  */


  process bwamem {

    tag "$replicate_id"
    publishDir "${params.aligned}", mode: 'copy'

    input:
    set val(replicate_id), file(reads) from normal_pairs_ch_2
    set val(replicate_id_2), file(reads_2) from tumor_pairs_ch_2
    file(bwa_index) from reference_index_ch_2

    output:
    set val(replicate_id), file("${replicate_id}.sorted.bam") into normal_aligned_ch
    set val(replicate_id_2), file("${replicate_id_2}.sorted.bam") into tumor_aligned_ch

    script:
    rg_n="\'@RG\\tID:${replicate_id}\\tSM:${params.sm}normal\\tPU:unit1\\tPL:illumina\\tLB:lib1\'"
    rg_t="\'@RG\\tID:${replicate_id_2}\\tSM:${params.sm}tumor\\tPU:unit1\\tPL:illumina\\tLB:lib1\'"
    """
    bwa mem -t 15 \
    -R $rg_n \
    $params.reference \
    $reads \
    | samtools sort -O bam - > ${replicate_id}.sorted.bam &

    bwa mem -t 15 \
    -R $rg_t \
    $params.reference \
    $reads_2 \
    | samtools sort -O bam - > ${replicate_id_2}.sorted.bam

    """

  }



  process creating_unmapped{

    tag "$replicate_id"
    publishDir "${params.aligned}", mode: 'copy'

    input:
    file reference from params.reference
    path refDir from params.referenceDir
    set val(replicate_id), file(read_fwd), file(read_rev) from separate_normal
    set val(replicate_id_2), file(read_fwd_2), file(read_rev_2) from separate_tumor

    output:
    tuple val(replicate_id), file("${replicate_id}.unaligned.bam") into normal_unmapped_ch
    tuple val(replicate_id_2), file("${replicate_id_2}.unaligned.bam") into tumor_unmapped_ch

    script:
    """
    picard FastqToSam \
    F1=$read_fwd \
    F2=$read_rev \
    R=$reference \
    O=${replicate_id}.unaligned.bam \
    SM=${params.sm}normal \
    RG=$replicate_id \
    PU=unit1 \
    PL=illumina \
    LB=lib1 &


    picard FastqToSam \
    F1=$read_fwd_2 \
    F2=$read_rev_2 \
    R=$reference \
    O=${replicate_id_2}.unaligned.bam \
    SM=${params.sm}tumor \
    RG=$replicate_id_2 \
    PU=unit1 \
    PL=illumina \
    LB=lib1
    """
  }

  process merging_exomes{

    publishDir "$params.outdir", mode:"copy"

    input:
    path index from reference_index_ch
    path dict from reference_dict_ch
    path reference from params.reference
    set val(replicate_id), path(aligned) from normal_aligned_ch
    set val(unaligned_id), path(unaligned) from normal_unmapped_ch
    set val(replicate_id_2), path(aligned_2) from tumor_aligned_ch
    set val(unaligned_id_2), path(unaligned_2) from tumor_unmapped_ch

    output:
    set val(replicate_id), file("${replicate_id}.merged.bam") into merged_normal_ch
    set val(replicate_id_2), file("${replicate_id_2}.merged.bam") into merged_tumor_ch

    script:
    """
    picard MergeBamAlignment \
    ALIGNED=$aligned \
    UNMAPPED=$unaligned \
    R=$reference \
    O=${replicate_id}.merged.bam &

    picard MergeBamAlignment \
    ALIGNED=$aligned_2 \
    UNMAPPED=$unaligned_2 \
    R=$reference \
    O=${replicate_id_2}.merged.bam
    """
  }

  process marking_duplicates_exomes{

    publishDir "$params.outdir", mode:"copy"

    input:
    set val(replicate_id), path(merged) from merged_normal_ch
    set val(replicate_id_2), path(merged_2) from merged_tumor_ch

    output:
    set val(replicate_id), file("${replicate_id}.merged.md.bam") into normal_duplicates_ch
    set val(replicate_id_2), file("${replicate_id_2}.merged.md.bam") into tumor_duplicates_ch
    file('*.marked_dup_metrics.txt')

    script:
    """
    picard MarkDuplicates \
    I=$merged \
    O=${replicate_id}.merged.md.bam \
    M=${replicate_id}.marked_dup_metrics.txt &

    picard MarkDuplicates \
    I=$merged_2 \
    O=${replicate_id_2}.merged.md.bam \
    M=${replicate_id_2}.marked_dup_metrics.txt

    """
  }


  process sorting_exomes{

    publishDir "$params.outdir", mode:"copy"

    input:
    set val(replicate_id), file(deduplicated) from normal_duplicates_ch
    set val(replicate_id_2), file(deduplicated_2) from tumor_duplicates_ch

    output:
    set val(replicate_id), file("${replicate_id}.sorted.RG.bam") into sorted_normal_ch
    set val(replicate_id_2), file("${replicate_id_2}.sorted.RG.bam") into sorted_tumor_ch

    script:
    """
    picard SortSam \
    I=$deduplicated \
    O=${replicate_id}.sorted.RG.bam \
    SORT_ORDER=coordinate &

    picard SortSam \
    I=$deduplicated_2 \
    O=${replicate_id_2}.sorted.RG.bam \
    SORT_ORDER=coordinate

    """
  }

  sorted_normal_ch.into { sorted_normal_ch_1; sorted_normal_ch_2 }
  sorted_tumor_ch.into { sorted_tumor_ch_1; sorted_tumor_ch_2 }

  process recalibration_exomes{

    publishDir "$params.outdir", mode:"copy"

    input:
    path index from reference_index_ch
    path dict from reference_dict_ch
    path reference from params.reference
    path data from params.dataDir
    set val(replicate_id), file(bam) from sorted_normal_ch_1
    set val(replicate_id_2), file(bam_2) from sorted_tumor_ch_1


    output:
    set val(replicate_id), file("${replicate_id}.recal_data.table") into recalibration_normal_ch
    set val(replicate_id_2), file("${replicate_id_2}.recal_data.table") into recalibration_tumor_ch


    script:
    """
    gatk BaseRecalibrator \
    -I $bam \
    -R $reference \
    --known-sites $params.wgs  \
    -O ${replicate_id}.recal_data.table \
    --disable-sequence-dictionary-validation &

    gatk BaseRecalibrator \
    -I $bam_2 \
    -R $reference \
    --known-sites $params.wgs \
    -O ${replicate_id_2}.recal_data.table \
    --disable-sequence-dictionary-validation
    """
  }


  process recalibration_BQSR_exomes{

    publishDir "$params.outdir", mode:"copy"

    input:
    path index from reference_index_ch
    path dict from reference_dict_ch
    path reference from params.reference
    set val(replicate_id), file(recalibrated) from recalibration_normal_ch
    set val(replicate_id_sort), file(bam) from sorted_normal_ch_2
    set val(replicate_id_2), file(recalibrated_2) from recalibration_tumor_ch
    set val(replicate_id_sort_2), file(bam_2) from sorted_tumor_ch_2


    output:
    set val(replicate_id), file("${replicate_id}.recalibration.bam"), file("${replicate_id}.recalibration.bai") into bqsr_normal_ch
    set val(replicate_id), file("${replicate_id_2}.recalibration.bam"), file("${replicate_id_2}.recalibration.bai") into bqsr_tumor_ch


    script:
    """
    gatk ApplyBQSR \
    -R $reference \
    -I $bam \
    --bqsr-recal-file $recalibrated \
    -O ${replicate_id}.recalibration.bam &

    gatk ApplyBQSR \
    -R $reference \
    -I $bam_2 \
    --bqsr-recal-file $recalibrated_2 \
    -O ${replicate_id_2}.recalibration.bam

    """
  }


 bqsr_normal_ch.into { bqsr_normal_ch_1; bqsr_normal_ch_2; bqsr_normal_ch_3 }
 bqsr_tumor_ch.into { bqsr_tumor_ch_1; bqsr_tumor_ch_2; bqsr_tumor_ch_3 }


 /*
 *  GENOTYPING, FILTRATION AND ANNOTATION
 */

  process genotyping{

    publishDir "$params.geno", mode:"copy"

    input:
    set val(replicate_id), path(normal), path(normal_index) from bqsr_normal_ch_1
    set val(replicate_id_2), path(tumor), path(tumor_index) from bqsr_tumor_ch_1
    path index from reference_index_ch
    path dict from reference_dict_ch
    path reference from params.reference

    output:
    path "*.{bcf.gz,bcf.gz.stats,bcf.gz.tbi}" into genotyping_ch

    script:
    id = "${replicate_id}_${replicate_id_2}"
    """
    samtools index $normal & samtools index $tumor

    gatk Mutect2 \
    -R $reference \
    -I $tumor \
    -I $normal \
    -normal ${params.sm}normal \
    -O ${id}.bcf.gz

    """
  }



process estimate_contamination{

  publishDir "$params.geno", mode:"copy"

  input:
  set val(replicate_id), path(normal), path(normal_index) from bqsr_normal_ch_2
  set val(replicate_id_2), path(tumor), path(tumor_index) from bqsr_tumor_ch_2
  path(data) from params.dataDir

  output:
  file("contamination.filt.table") into contamination_ch


  script:

  """
  gatk GetPileupSummaries --disable-sequence-dictionary-validation \
  -I $tumor \
  -V ${data}/ALL.wgs.biallelic.exons.af.vcf.gz \
  -L ${data}/ALL.wgs.biallelic.exons.af.vcf.gz \
  -O ${params.geno}/pileups.filt.table &

  gatk GetPileupSummaries --disable-sequence-dictionary-validation \
  -I $normal \
  -V ${data}/ALL.wgs.biallelic.exons.af.vcf.gz \
  -L ${data}/ALL.wgs.biallelic.exons.af.vcf.gz \
  -O ${params.geno}/pileups.normal.filt.table &
  wait

  gatk CalculateContamination \
  -I ${params.geno}/pileups.filt.table \
  -matched ${params.geno}/pileups.normal.filt.table \
  --tumor-segmentation $params.segments \
  -O contamination.filt.table

  """
}


process orientation_bias{

  publishDir "$params.geno", mode:"copy"

  input:
  set val(replicate_id), path(tumor), path(tumor_index) from bqsr_tumor_ch_3
  path index from reference_index_ch
  path dict from reference_dict_ch
  path reference from params.reference

  output:
  file("tumor.artifacts.tar.gz") into orientation_bias_ch

  script:
  """
  gatk CollectF1R2Counts \
  -I $tumor \
  -O tumor.F1R2.tar.gz \
  -R $reference

  gatk LearnReadOrientationModel \
  -I tumor.F1R2.tar.gz \
  -O tumor.artifacts.tar.gz

  """
}


  process filter_calls{

    publishDir "$params.geno", mode:"copy"

    input:
    path index from reference_index_ch
    path dict from reference_dict_ch
    path reference from params.reference
    set path(bcf), path(index), path(stats) from genotyping_ch
    path(contamination) from contamination_ch
    path(artifacts) from orientation_bias_ch


    output:
    path("*.{vcf.gz,tsv,vcf.gz.tbi}") into filtered_ch

    script:
    """
    gatk FilterMutectCalls \
    -R $reference \
    -V $bcf \
    --contamination-table $contamination \
    --tumor-segmentation $params.segments \
    --orientation-bias-artifact-priors $artifacts \
    -O ${bcf.simpleName}.filtered.vcf.gz


    bcftools filter -i'FILTER="PASS"' -O z ${bcf.simpleName}.filtered.vcf.gz

    gunzip ${bcf.simpleName}.filtered.vcf.gz
    """
  }



  process annotate_calls{

  publishDir "$params.geno", mode:"copy"

  input:
  set path(variants), path(stats), path(index) from filtered_ch
  path index from reference_index_ch
  path dict from reference_dict_ch
  path reference from params.reference

  output:
  file("*.{vcf.gz, tbi, idx}") into annotated_ch

  script:
  """
  gatk Funcotator \
  --variant $variants \
  --reference $reference \
  --ref-version hg38 \
  --data-sources-path $params.funcotator_data \
  --output ${variants.simpleName}.funcotated.vcf \
  --output-file-format VCF

  bgzip -c ${variants.simpleName}.funcotated.vcf > ${variants.simpleName}.funcotated.vcf.gz

  tabix -p vcf ${variants.simpleName}.funcotated.vcf.gz

  """
  }
