process gatk4_createsequencedict {

    input:
    path fasta

    output:
    path "*.dict" , emit: dict

    script:
    """
    gatk CreateSequenceDictionary -REFERENCE $fasta -OUTPUT ${fasta.baseName}.dict
    """
}

process gatk4_createintervallist {
  input:
    path(reffaidx)
    path(refdict)

  output:
    path("*.interval_list")


  script:
    def outfile = "${refdict.baseName}.interval_list"

  """
  cat $reffaidx | awk '{OFS="\t";print \$1,0,\$2}' > ref.bed
  gatk BedToIntervalList -I ref.bed -O $outfile -SD $refdict
  """
}

process gatk_scatterintervals {
  input:
    path(intervallist)
    val(chunksize)

  output:
    path("scatter_temp*.interval_list")

  script:

  """
  mkdir -p scatter
  gatk IntervalListTools -I $intervallist --SCATTER_CONTENT $chunksize -O scatter
  for f in \$(find scatter -name '*.interval_list');do 
    tn=\$(echo \$f | sed 's:/:_:g');
    cp \$f \$tn
  done
  """
}


process gatk_baserecalibrator {
  input:
    tuple val(sample), path(bam), path(bai), path(interval)
    path(genome)
    path(index)
    path(dict)
    tuple val(snpmeta), path(snps), path(snpsi)
    tuple val(indelmeta), path(indels), path(indelsi)    


  output:
    tuple val(sample), val(interval.baseName), path("*.recal"), emit: recal

  script:
    def args = task.ext.args ?: ''
    def outfile = "${bam.baseName}.${interval.baseName}.recal"

  """
  gatk --java-options "-Xmx${task.memory.giga}G" \
      BaseRecalibrator \
      -R $genome \
      -I $bam \
      --use-original-qualities \
      -O $outfile \
      --known-sites $snps \
      --known-sites $indels \
      -L $interval
  """
}


process gatk_gatherbqsrreports {

  publishDir "$params.outdir/bqsr", mode: 'copy'

  input:
    tuple val(sample), path(reportlist)
    val(round)

  output:
    tuple val(sample), path("*merged.recal")

  script:
    def outfile = "${sample}.${round}.merged.recal"
    def infiles = reportlist.join(" -I ")

  """
  gatk --java-options "-Xmx${task.memory.giga}G" \
     GatherBQSRReports \
      -I ${infiles} \
      -O ${outfile}
  """
}

// Note that parameters for this are a bit unintuitive. See https://gatk.broadinstitute.org/hc/en-us/community/posts/360056320991-options-used-during-base-recalibration
// Apparently it is fine to quantize quals into these crude bins 10,20,30
// --use-original-qualities allows us to apply BQSR in a subsequent round
//
process gatk_applybqsr {
  input:
    tuple val(sample), path(bam), path(bai), path(interval), path(recalreport)
    path(genome)
    path(index)
    path(dict)

  output:
    tuple val(sample), path("*recal.bam")

  script:
    def outfile = "${bam.baseName}.${interval.baseName}.recal.bam"

  """
  gatk --java-options "-Xms${task.memory.giga}G" \
      ApplyBQSR \
      -R ${genome} \
      -I ${bam} \
      -O ${outfile} \
      -L ${interval} \
      -bqsr ${recalreport} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
  """

}

process gatk_mergebams {

  publishDir "$params.outdir/bqsr", mode: 'copy'

  input:
    tuple val(sample), path(bamlist)

  output:
    tuple val(sample), path("*recal.merged.bam"), path("*recal.merged.*bai")

  script:
    def outfile = "${sample}.recal.merged.bam"
    def infiles = bamlist.join(" --INPUT ")

  """
  gatk --java-options "-Dsamjdk.compression_level=5 -Xms${task.memory.giga - 1}G" \
      GatherBamFiles \
      --INPUT ${infiles} \
      --OUTPUT ${outfile} \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true
  """
}

process gatk_analyzecovariates {

  publishDir "$params.outdir/bqsr", mode: 'copy'

  input:
    tuple val(sample), path(report)

  output:
    tuple val(sample), path("*.pdf")

  script:
    def outfile = "${sample}.recal.report.pdf"

  """
  gatk --java-options "-Xmx${task.memory.giga - 1}G" \
    AnalyzeCovariates \
    -bqsr ${report} \
    -plots ${outfile}
  """

}

