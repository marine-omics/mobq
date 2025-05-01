nextflow.enable.dsl=2

include { gatk4_createsequencedict; gatk4_createintervallist; gatk_scatterintervals;  gatk_baserecalibrator; gatk_gatherbqsrreports; gatk_applybqsr; gatk_mergebams; gatk_analyzecovariates } from './modules/gatk.nf'
include { sidx; faidx } from './modules/samtools.nf'
include { tabix } from './modules/bcftools.nf'

// Command line parameters
// params.snps = null
// params.indels = null
params.bams = null



workflow do_bqsr {
  take:
    indexed_bams
    genome_fasta
    genome_dict
    ch_snps
    ch_indels

  main:

    ch_bamcollection = indexed_bams.map{m,b,i -> b} | collect
    ch_baicollection = indexed_bams.map{m,b,i -> i} | collect

    genome_fai = faidx(genome_fasta) | collect

    genome_intervals = gatk4_createintervallist(genome_fai,genome_dict)
    ch_gatk_scatter_intervals = gatk_scatterintervals(genome_intervals,params.gatk_chunksize) | flatten

    ch_baserecalibrator_inputs = indexed_bams.combine(ch_gatk_scatter_intervals)


    ch_recal_tables = gatk_baserecalibrator(ch_baserecalibrator_inputs,genome_fasta,genome_fai,genome_dict,ch_snps,ch_indels)

    ch_merged_reports = ch_recal_tables.groupTuple().map{ s,interval,recal -> [s,recal]} | gatk_gatherbqsrreports

  //  ch_baserecalibrator_inputs.view()

//    ch_merged_reports.view()

    //ch_merged_reports.join(ch_baserecalibrator_inputs,remainder:true).view()
    ch_applybqsr_inputs = ch_baserecalibrator_inputs.combine(ch_merged_reports, by: 0)

    ch_recal_bams = gatk_applybqsr(ch_applybqsr_inputs,genome_fasta,genome_fai,genome_dict)

    ch_final_bams = ch_recal_bams.groupTuple() | gatk_mergebams

    ch_merged_reports | gatk_analyzecovariates

}


workflow {
  genome_fasta = Channel.fromPath(file(params.genome, checkIfExists:true)) | collect
  genome_dict = gatk4_createsequencedict(genome_fasta) | collect

  snps_vcf = ["snps",file(params.snps, checkIfExists:true)]
  indels_vcf = ["indels",file(params.indels, checkIfExists:true)]

  ch_mapped_marked_bams = extract_samples_bam(file(params.bams,checkIfExists:true))

  ch_mapped_marked_bais = ch_mapped_marked_bams | sidx

  ch_bbai = ch_mapped_marked_bams.join(ch_mapped_marked_bais)

  ch_vcf = Channel.fromList( [snps_vcf, indels_vcf])
  tbi_l = ch_vcf | tabix 

//  ch_vcfi = ch_vcf.join(ch_tbi)
  ch_vcfi = ch_vcf.join(tbi_l).branch {
    m,v,i -> 
    snps: m=="snps"
    indels: m=="indels"
  }.
  set { result }


  ch_snps = result.snps.collect()
  ch_indels = result.indels.collect()

  do_bqsr(ch_bbai,genome_fasta,genome_dict,ch_snps,ch_indels)

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def resolve_path(pathstring){
  if(pathstring =~ /^\//){
    pathstring
  } else {
    "${params.base_path}/${pathstring}"
  }
}

def extract_samples_bam(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
    .map{ row -> 
      def sample = row.sample
      def bam     = file(resolve_path(row.bam), checkIfExists: true)

      [sample,bam]
    }
}
