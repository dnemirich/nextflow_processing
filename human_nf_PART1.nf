
// INPUT FILES


params.filtered = "$baseDir/input/sample1/transcriptomes/filtered_transcriptomes/*{1,2}_filtered.noadapt.fastq.gz"
params.qc = "$baseDir/qc/sample1"
params.reference = "$baseDir/reference/hg38.analysisSet.fa"
params.referenceDir = "$baseDir/reference"
params.dataDir = "$baseDir/data"
params.sortDir = "$baseDir/sorted/sample1"
params.sm = "sample1_tumor_rna"
params.merged = "$baseDir/merged/sample1"
params.genotyping = "$baseDir/genotyping/sample1"



// OUTPUT CONFIGURATION

params.outdir = './results'



 Channel
 .fromFilePairs(params.filtered, flat: true)
 .set { separate_reads }

 Channel
 .fromFilePairs( params.filtered )
 .ifEmpty { error "Cannot find any reads matching: ${params.filtered}" }
 .into { read_pairs_ch; read_pairs_ch2 }




 /*
  * QC STEP FOR TRIMMED READS
  */

  process trimmed_fastqc{
    publishDir "${params.qc}/trimmed", mode: 'copy'

    input:
    set val(prefix), file(reads) from read_pairs_ch

    output:
    file "*.{zip,html,txt}" into fastqc_trimmed_results

    script:

    """
    fastqc -t 10 $reads
    """

  }


 /*
  *  REFERENCE PREPARATION STEP
  */


  process reference_preparation_1{
    
    tag "$params.reference"
    publishDir path: params.referenceDir, mode: 'copy'

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
    publishDir path: params.referenceDir, mode: 'copy'

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
    publishDir path: params.referenceDir, mode: 'copy'
    
    input:
    path reference from params.reference

    output:
    file "*.{amb,ann,bwt,pac,sa}" into reference_index_ch_2

    script:
    """
    bwa index $reference
    """
  }


  process reference_preparation_4{

    tag "$params.reference"

    input:
    path referenceDir from params.referenceDir
    path reference from params.reference

    output:

    file "star" into star_reference_ch

    script:
    """
    mkdir star

    STAR --runThreadN 20 --runMode genomeGenerate \
    --genomeDir star/ \
    --genomeFastaFiles $reference \
    --sjdbGTFfile ${params.dataDir}/gencode.v31.chr_patch_hapl_scaff.annotation.gtf \
    --sjdbOverhang 99
    """
  }

  
 /*
 *  TRANSCRIPTOMES PROCESSING
 */


 process transcriptomes_STAR_alignment{

   tag "$replicate_id"

   input:
   path star from star_reference_ch
   path reference from params.reference
   tuple  val(replicate_id), path(reads) from read_pairs_ch2

   output:
   tuple val(replicate_id), path('*_Aligned.out.bam') into aligned_bam_ch

   script:
   """
   STAR --genomeDir $star \
        --readFilesIn $reads \
        --runThreadN 20 \
        --readFilesCommand zcat \
        --outFileNamePrefix ${replicate_id}_ \
        --outSAMtype BAM Unsorted  \
   """
 }


 process sorting_mapped{

   tag "$replicate_id"
   publishDir "$params.sortDir", mode:"copy"

   input:
   tuple val(replicate_id), path(bam) from aligned_bam_ch

   output:
   tuple val(replicate_id), path('*.sorted.bam') into sorted_bam_ch

   script:
   """
   picard SortSam I=$bam \
   O=${replicate_id}.sorted.bam \
   SO=queryname
   """
 }


 process merging_unmapped{

   tag "$replicate_id"

   input:
   path reference from params.reference
   path refDir from params.referenceDir
   set val(replicate_id), file(read_fwd), file(read_rev) from separate_reads

   output:
   tuple val(replicate_id), path('*.unaligned.bam') into merged_bam_ch

   script:
   """
   picard FastqToSam \
   F1=$read_fwd \
   F2=$read_rev \
   R=$reference \
   O=${replicate_id}.unaligned.bam \
   SM=${params.sm} \
   RG=$replicate_id \
   LB=lib1 \
   PL=illumina \
   PU=unit1
   """
 }


 process sorting_unmapped{

   tag "$replicate_id"
   publishDir "$params.sortDir", mode:"copy"

   input:
   tuple val(replicate_id), path(bam) from merged_bam_ch

   output:
   tuple val(replicate_id), path('*.sorted.bam') into sorted_unmapped_bam_ch
   val true into done_ch

   script:
   """
   picard SortSam \
   I=$bam \
   O=${replicate_id}.unaligned.sorted.bam \
   SO=queryname
   """
 }

 process merging_mapped_n_unmapped{
   
   input:
   val flag from done_ch

   """
   bash $baseDir/merging.sh -w $params.sortDir -r $params.reference -o $params.merged
   """
 }


 Channel
  .fromPath("${params.merged}/*.merged.bam")
  .set { merged_mapped_n_unmapped_ch }


 process marking_duplicates{
   
   publishDir "$params.outdir", mode:"copy"

   input:
   path(merged) from merged_mapped_n_unmapped_ch

   output:
   path('*.md.bam') into duplicates_ch
   file('*.marked_dup_metrics.txt')

   script:
   """
   picard MarkDuplicates \
   I=$merged \
   O=${merged.simpleName}.md.bam \
   M=${merged.simpleName}.marked_dup_metrics.txt
   """
 }


 process splitting_reads{

   input:
   path index from reference_index_ch
   path dict from reference_dict_ch
   path reference from params.reference
   path(duplicates) from duplicates_ch

   output:
   set path("*.bai"), path('*.split.bam') into splitted_ch

   script:
   """
   gatk SplitNCigarReads \
   -R $reference \
   -I $duplicates \
   -O ${duplicates.simpleName}.split.bam
   """
 }


 splitted_ch.into { splitted_ch_1; splitted_ch_2 }


 process recalibrator{
   
   publishDir "$params.outdir", mode:"copy"
   
   input:
   path index from reference_index_ch
   path dict from reference_dict_ch
   path reference from params.reference
   set path(bai), path(bam) from splitted_ch_1
   path data from params.dataDir

   output:
   path("*.recal_data.table") into recalibration_ch

   script:
   """
   gatk BaseRecalibrator \
   -I $bam \
   -R $reference \
   --known-sites $data/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf \
   -O ${bam.simpleName}.recal_data.table \
   --disable-sequence-dictionary-validation
   """
 }


 process recalibration_BQSR{
   
   publishDir "$params.outdir", mode:"copy"
   
   input:
   path index from reference_index_ch
   path dict from reference_dict_ch
   path reference from params.reference
   path(recalibrated) from recalibration_ch
   set path(bai), path(bam) from splitted_ch_2

   output:
   path("*.recal.bam") into bqsr_ch

   script:
   """
   gatk ApplyBQSR \
   -R $reference \
   -I $bam\
   --bqsr-recal-file $recalibrated \
   -O ${recalibrated.simpleName}.recal.bam

   """
 }


process haplotype_caller{

  publishDir "$params.genotyping", mode:"copy"

  input:
  path reference from params.reference
  path index from reference_index_ch
  path dict from reference_dict_ch
  path bam from bqsr_ch
  path geno from params.genotyping

  output:
  set file("*.vcf"), file("*.vcf.idx") into haplotypecaller_ch


  script:
  """
  gatk HaplotypeCaller \
  -R $reference \
  -I $bam \
  -O ${bam.simpleName}.vcf
  """
}
