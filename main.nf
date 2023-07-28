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
  bbmerge-auto.sh in1=${read_1} in2=${read_2} out=bbmerged.fastq.gz outu1="unbbmerged_${read_1}" outu2="unbbmerged_${read_2}" &> bbmerge.log
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

process MOTUS_PROFILE {
  input:
    tuple path(read_1), path(read_2)
    val(meta_sample_name)
  output:
  
  script:
  """
  motus profile -f ${read_1} -r ${read_2} -t ${task.cpus}
  """
}

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

        // Assembly for determining insert size
        //   If data had higher coverage and resources allowed, assembly based taxonomy assessment could be performed
        //  Testing with mmseq2 Swissprot database on assembly yielded 'Bacteroides_thetaiotaomicron' but low confidence w/ alignment to 48 reference genome sample for species yielding low alignment rates (<6%)
        //  Resource limitation: mainly reference database on disk limits
        MEGAHIT( BBDUK.out.trimmed_reads, params.meta_sample_name )

        MOTUS_SETUP() // Run Motus setup once to motus database

        MOTUS_PROFILE( BBDUK.out.trimmed_reads, params.meta_sample_name )

        MOTUS_PROFILE.out.top_species | GET_NCBI_TAXONOMY
        
        BUILD_BOWTIE2_REFERENCE( 
          GET_NCBI_TAXONOMY.out.taxonomy, 
          params.number_genome_representatives 
          )

        ALIGN_BOWTIE2( 
          BUILD_BOWTIE2_REFERENCE.out.reference,
          BBDUK.out.trimmed_reads
          )

        ALIGN_BOWTIE2.out.sam | SAMTOOLS_SORT | SAMTOOLS_INDEX | EXTRACT_PER_BASE_QUALITY

        GENERATE_REPORTS(
          MEGAHIT.out.log, // For insert size
          TRIMMED_FASTQC.out.report_zip, // For plot, per position quality
          EXTRACT_PER_BASE_QUALITY.out.per_base_position_quality, // For plot, per postion alignment quality
          GET_NCBI_TAXONOMY.out.taxonomy // For text report, species name and ncbi taxid 
          BUILD_BOWTIE2_REFERENCE.out.full_reference_list // For text report, species genome accessions
          )
        
}