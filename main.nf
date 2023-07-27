include { FASTQC as RAW_FASTQC } from "./modules/FASTQC.nf"
include { FASTQC as TRIM_FASTQC} from "./modules/FASTQC.nf"

process BBDUK {
  publishDir 'bbduk_out'

  input:
    tuple path(read_1), path(read_2)
    path(adapters_ref)
  output:
    tuple path("${read_1.simpleName}-trimmed.fastq.gz"), path("${read_2.simpleName}-trimmed.fastq.gz"), emit: trimmed_reads
  
  script:
  """
  bbduk.sh in=${read_1} in2=${read_2} out1=${read_1.simpleName}-trimmed.fastq.gz \
         out2=${read_2.simpleName}-trimmed.fastq.gz ref=${adapters_ref} ktrim=l k=17 ftm=5 qtrim=rl \
         trimq=10 mlf=0.5 maxns=0
  """
}

process BBMERGE {
    publishDir "bbmerge"

  input:
      tuple path(read_1), path(read_2)

  output:
    path("bbmerged.fastq.gz"), emit: bbmerged_reads
    path("unbbmerged_${read_1}"), emit: unbbmerged_1
    path("unbbmerged_${read_2}"), emit: unbbmerged_2
    path("bbmerge.log"), emit: bbmerge_log
  
  script:
  """
  bbmerge.sh in1=${read_1} in2=${read_2} out=bbmerged.fastq.gz outu1="unbbmerged_${read_1}" outu2="unbbmerged_${read_2}" > bbmerge.log
  """
}

process MULTIQC {
  publishDir "MultiQC"

  input:
    path multiqc_files, stageAs: "?/*"
  
  output:
    path("multiqc_report.html"), emit: multiqc_report
  
  script:
  """
  multiqc .
  """
}

process MEGAHIT {
    publishDir "megahit_out"

  input:
    tuple path(read_1), path(read_2)
    val(meta_sample_name)
  output:
    path("${meta_sample_name}"), emit: contig_dir
  
  script:
  """
  megahit -1 ${read_1} -2 ${read_2} \
        -o ${meta_sample_name} --min-contig-len 500
  """
}
/*

process GET_CAT_DB {
    storeDir "CAT_DB"

  input:
    val(target_db_name)
    
  output:
    path(target_db_name)

  script:
  """
  wget 'https://tbb.bio.uu.nl/bastiaan/CAT_prepare/${target_db_name}.tar.gz'
  tar -xvzf ${target_db_name}.tar.gz
  """
}

process CAT {
  input:
    path(contig_fa)
    path(CAT_DB)
  
  output:
  
  script:
  """
  CAT contigs -c ${contig_fa} -d ${CAT_DB} \
            -t CAT_prepare_20200618/2020-06-18_taxonomy/ -p sample-1-genes.faa \
            -o sample-1-tax-out.tmp -n 15 -r 3 --top 4 --I_know_what_Im_doing
  """
}

process RENAME_AND_SUMMARIZE_CONTIGS {
  input:
    path(contig_dir)
  output:
    path("${meta_sample_name}.fasta"), emit: renamed_contigs
    path("assembly-summaries.tsv"), emit: assembly_sumary
  
  script:
  """
  bit-rename-fasta-headers -i ${contig_dir}/final.contigs.fa -w c_${meta_sample_name} -o ${meta_sample_name}.fasta
  bit-summarize-assembly -o assembly-summaries.tsv ${meta_sample_name}.fasta

  """
}
*/
workflow {
    main:
        ch_reads = Channel.fromPath(params.reads_dir + "/*.fastq.gz") | collect
        
        // Database setups, only need to run once as results with be stored for future runs
        // GET_CAT_DB(params.cat_db_name)

        // Actual pipeline
        RAW_FASTQC(ch_reads, "raw")
        BBDUK(ch_reads, params.bbduk_adapters)
        TRIM_FASTQC(BBDUK.out.trimmed_reads, "trimmed")

        MULTIQC(RAW_FASTQC.out.fastqc_zip | mix(TRIM_FASTQC.out.fastqc_zip) | collect)

        MEGAHIT( BBDUK.out.trimmed_reads, params.meta_sample_name)
        BBMERGE( BBDUK.out.trimmed_reads)

        // CAT( RENAME_AND_SUMMARIZE_CONTIGS.out.renamed_contigs, GET_CAT_DB.out.CAT_DB )
}