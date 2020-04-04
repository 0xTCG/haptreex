// 786

params.bam = ""
params.repo = "."
params.data = "${params.repo}/broad-grch37"
params.tools = "${params.repo}/tools"
params.output = "output"
params.overhang = 75
params.sample = "sample"

log.info """\
========================================================================
=                          GATK RNA-seq caller                         =
========================================================================
"""

gatk = "${params.tools}/gatk4/gatk"
gatk_extra = "--java-options \"${params.java} -Xms6144m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10\""
picard = "java ${params.java} -jar ${params.tools}/picard.jar"

input_bam_ch = Channel.fromPath(params.bam).map { file -> tuple(file.baseName, file) }

process 'RevertSam' {
  input:  set id, file(bam) from input_bam_ch
  output: set id, file('reverted.bam') into revert_ch, revert_merge_ch
  script:
    """
    ${gatk} RevertSam \
      --INPUT ${bam} \
      --OUTPUT reverted.bam \
      --VALIDATION_STRINGENCY SILENT \
      --ATTRIBUTE_TO_CLEAR FT \
      --ATTRIBUTE_TO_CLEAR CO \
      --SORT_ORDER queryname
    """
}

process 'SamToFastq' {
  input:  set id, file(bam) from revert_ch
  output: set id, file("sample_1.fq.gz"), file("sample_2.fq.gz") into fq_ch
  script:
    """ 
    ${gatk} SamToFastq \
      --INPUT ${bam} \
      --VALIDATION_STRINGENCY SILENT \
      --FASTQ sample_1.fq.gz \
      --SECOND_END_FASTQ sample_2.fq.gz
    """
}

process 'MapFastq' {
  tag "$id"
  
  module 'star'
  cpus 16
  memory { 40.GB * task.attempt }

  input:  set id, file(fq1), file(fq2) from fq_ch
  output: set id, file("star_map.Aligned.sortedByCoord.out.bam") into star_ch
  script:
    """
    STAR \
      --genomeDir ${params.data}/STAR2_${params.overhang} \
      --runThreadN 16 \
      --readFilesIn ${fq1} ${fq2} \
      --readFilesCommand "gunzip -c" \
      --sjdbOverhang ${params.overhang} \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --outSAMattrRGline ID:id LB:rnaseq PL:illumina SM:${id} \
      --outFileNamePrefix star_map.
    """
}

process 'MarkDuplicates' {
  tag "$id"
  echo true
  publishDir "${params.output}/$id"

  module 'sambamba'
  cpus 8
  time 24.hours
  memory { 48.GB * task.attempt }

  input:  set id, file (bam) from star_ch
  output: set id, file ('dedupped.bam') into dedup_ch
  script:
    """
    sambamba markdup \
      -t ${task.cpus} \
      -p \
      ${bam} \
      dedupped.bam
    """
}

process 'SplitNCigar' {
  tag "${id}.${interval}"
  echo true

  input:
    set id, file(bam) from dedup_ch
    each interval from Channel.from(1..25)
  output:
    set id, interval, file("split.${interval}.bam") into split_part_ch
  script:
    """
    ${gatk} SplitNCigarReads \
      -I ${bam} \
      -O split.${interval}.bam \
      -R ${params.data}/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
      -L ${params.data}/star.gencode.v19.transcripts.patched_contigs.exons.interval_list.scattered/${interval}.part.interval_list
    """
}

splits = split_part_ch.map { id, intvl, f -> [id, [intvl, f]] }.groupTuple (size: 25)
process 'BaseRecalibrator' {
  tag "$id"
  echo true

  input:  set id, files from splits
  output: set id, file('recal.out'), files into recal_ch

  script:
    """
    ${gatk} ${gatk_extra} BaseRecalibrator \
      -R ${params.data}/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
      -O recal.out \
      --use-original-qualities \
      -known-sites ${params.data}/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf \
      -known-sites ${params.data}/Mills_and_1000G_gold_standard.indels.b37.sites.vcf \
      -known-sites ${params.data}/Homo_sapiens_assembly19_1000genomes_decoy.known_indels.vcf \
      ${files.collect { _, j -> "-I $j" }.join(' ') }
    """
}

recals = recal_ch.flatMap { id, recal, fls -> fls.collect { intvl, f -> [id, intvl, f, recal] } }
process 'ApplyBQSR' {
  tag "${id}.${interval}"
  echo true
  publishDir "${params.output}/$id/bqsr"

  input:  set id, interval, file(bam), file(recal) from recals
  output: set id, interval, file("bqsr.${interval}.bam") into bqsr_ch, bcf_ch
  script:
    """
    ${gatk} ApplyBQSR \
      -R ${params.data}/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
      -I ${bam} \
      -O bqsr.${interval}.bam \
      --bqsr-recal-file ${recal} \
      --add-output-sam-program-record \
      --use-original-qualities
    """
}

process 'HaplotypeCaller' {
  tag "${id}.${interval}"
  echo true
  publishDir "${params.output}/${id}/hc"
  
  cpus 8
  
  input:  
    set id, interval, file(bam) from bqsr_ch
  output: 
    set id, interval, file("output.${interval}.vcf.gz"), file("output.${interval}.vcf.gz.tbi") into call_ch
  script:
    """
    ${gatk} ${gatk_extra} HaplotypeCaller \
      -R ${params.data}/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
      -I ${bam} \
      --dbsnp ${params.data}/Homo_sapiens_assembly19_1000genomes_decoy.dbsnp138.vcf \
      -O output.${interval}.vcf.gz \
      -dont-use-soft-clipped-bases \
      --native-pair-hmm-threads ${task.cpus} \
      --standard-min-confidence-threshold-for-calling 20 \
      --L ${params.data}/star.gencode.v19.transcripts.patched_contigs.exons.interval_list.scattered/${interval}.part.interval_list
    """
}

calls = call_ch.map { id, intvl, f, fidx -> [id, [f, fidx]] }.groupTuple (size: 25)
process 'VariantFiltration' {
  tag "$id"
  echo true
  publishDir "${params.output}/$id"

  input:  set id, files from calls
  output: set file('final.vcf'), file('final.vcf.idx')
  script:
    """
    ${gatk} MergeVcfs \
      ${files.collect { f, _ -> "-I $f" }.join(' ') } \
      -O merged.vcf

    ${gatk} VariantFiltration \
      --R ${params.data}/Homo_sapiens_assembly19_1000genomes_decoy.fasta \
      --V merged.vcf \
      -O final.vcf \
      --window 35 \
      --cluster 3 \
      --filter-name "FS" \
      --filter "FS > 30.0" \
      --filter-name "QD" \
      --filter "QD < 2.0"
    """
}

