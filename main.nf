include { FASTQC as RAW_FASTQC } from "./modules/FASTQC.nf"
include { FASTQC as TRIM_FASTQC} from "./modules/FASTQC.nf"

process BBDUK {
  publishDir "${params.output}/bbduk_out"

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
    publishDir "${params.output}/bbmerge"

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
  publishDir "${params.output}/MultiQC"

  input:
    path multiqc_files, stageAs: "?/*"
  
  output:
    path("multiqc_report.html"), emit: multiqc_report
    path("multiqc_data"), emit: data
  
  script:
  """
  multiqc .
  """
}

process MEGAHIT {
    publishDir "${params.output}/megahit_out"

  input:
    tuple path(read_1), path(read_2)
    val(meta_sample_name)
  output:
    path("${meta_sample_name}"), emit: contig_dir
    path("${meta_sample_name}/log"), emit: log
  
  script:
  """
  megahit -1 ${read_1} -2 ${read_2} \
        -o ${meta_sample_name} --min-contig-len 500
  """
}

// TODO: Replace downloadDB call with preloaded database in image; reduces first run setup time and locks in version of mOTU for reproducibilit
process MOTUS_PROFILE {
  publishDir "${params.output}/motus_out"

  input:
    tuple path(read_1), path(read_2)
    val(meta_sample_name)
  output:
    path("taxonomy_profile.txt"), emit: full_profile
    path("top_species.txt"), emit: top_species

  stub:
  """
  # Load precomputed results
  cp ${params.stub_motus_results} .

  # Create top species file
  extract_top_species.py
  """
  
  script:
  """
  motus downloadDB
  motus profile -f ${read_1} -r ${read_2} -t ${task.cpus} > taxonomy_profile.txt

  # Create top species file
  extract_top_species.py
  """
}

process FETCH_NCBI_GENOMES {
  storeDir "ncbi_genomes"

  input:
    val(top_species)
    val(number_genome_representatives)
  
  output:
    path("${top_species.replace(' ','_')}"), emit: genomes
    path("${top_species.replace(' ','_')}/accession_ids.txt"), emit: full_genomes_list
  
  script:
  species_dir = "${top_species.replace(' ','_')}"
  """
  mkdir ${species_dir}
  cd ${species_dir}
  get_genomes.sh '${top_species}' ${number_genome_representatives}
  """
}

process BUILD_BOWTIE2_REFERENCE {
  publishDir "${params.output}/bowtie2_references"

  input:
    val(top_species)
    path(genome), stageAs: "genomes"
  
  output:
    path("${species_dir}"), emit: reference
  
  script:
  species_dir = "${top_species.replace(' ','_')}"
  """
  build_bowtie.sh ${species_dir}
  mkdir ${species_dir}
  mv *.bt2 ${species_dir}
  """
}

process ALIGN_BOWTIE2 {
  publishDir "${params.output}/bowtie2_alignments"

  input:
    tuple path(read_1), path(read_2)
    path(reference)
    val(meta_sample_name)
  
  output:
    path("${meta_sample_name}_${reference}.sam"), emit: alignment
    path("${meta_sample_name}_alignment.log"), emit: log
  
  script:
  """
  bowtie2 -x ${reference}/${reference} \
      -1 ${read_1} \
      -2 ${read_2} \
      --no-unal \
      -p ${task.cpus} \
      -S ${meta_sample_name}_${reference}.sam &> ${meta_sample_name}_alignment.log
  """
}

process EXTRACT_PER_BASE_ALIGNMENT_SCORES {
  publishDir "${params.output}/quality_metrics"

  input:
    path(alignment) // sam format, unsorted
    val(sample_count_alignments) // How many alignment to generate the table based on, -1 denotes all
  
  output:
    path("per_positon_alignment_scores.csv"), emit: table
  
  script:
  """
  alignment_scores.py ${alignment} ${sample_count_alignments}
  """
}

process EXTRACT_PER_BASE_PHRED {
  publishDir "${params.output}/quality_metrics"

  input:
    tuple path(read_1), path(read_2)
  
  output:
    path("per_positon_quality_scores.csv"), emit: table
  
  script:
  """
  quality_scores.py ${read_1},${read_2}
  """
}

process GENERATE_REPORTS {
  publishDir "${params.output}/full_report"

  input:
    path(assembly_log), stageAs: 'assembly_log' // MEGAHIT.out.log, // For insert size mean and s.d.
    path(multiqc_data), stageAs: 'multiqc_data'
    path(fastqc_report)
    path(alignment_quality_table) // EXTRACT_PER_BASE_QUALITY.out.per_base_position_quality, // For plot, per postion alignment quality
    path(phred_quality_table) // EXTRACT_PER_BASE_QUALITY.out.per_base_position_quality, // For plot, per postion alignment quality
    path(full_reference_list) // BUILD_BOWTIE2_REFERENCE.out.full_reference_list // For text report, species genome accessions
    val(top_species)
    tuple path(read_1), path(read_2)
          
  output:
    path("report.txt"), emit: text_report
    path("Quality_scores.png"), emit: plots
  
  script:
  """

  # Extract read counts
  extract_total_sequences.py ${multiqc_data}/multiqc_general_stats.txt >> report.txt

  # Extract R1 Length, R1 Index, R2 Length, R2 Index
  echo R1 Length: \$(zcat ${read_1} | head -n 2 | tail -n 1 | tr -d '\n' | wc -c) >> report.txt
  echo R2 Length: \$(zcat ${read_2} | head -n 2 | tail -n 1 | tr -d '\n' | wc -c) >> report.txt
  echo R1 Index: \$(zcat ${read_1} | awk -F: '{print \$NF}' | head -n 1) >> report.txt
  echo R2 Index: \$(zcat ${read_2} | awk -F: '{print \$NF}' | head -n 1) >> report.txt
  
  # Extract insert size mean and s.d.
  echo "Insert Size Mean and S.D." >> report.txt
  extract_insert_size_stats.sh ${assembly_log} >> report.txt

  # Extract Nextera or TruSeq adapter content
  unzip *_fastqc.zip 
  sed -n '/>>Adapter Content/,/>>END_MODULE/p' */fastqc_data.txt > adapter_content.txt
  extract_adapter_content.py adapter_content.txt >> report.txt

  # Extract species idenfication information
  echo Species Idenfification using mOTUs: >> report.txt
  echo "  Top Species: " ${top_species} >> report.txt
  echo "  Reference Genomes For Top Species:" >> report.txt
  awk '{print "    - " \$0}' ${full_reference_list} >> report.txt

  # Generate Plots
  ## Per Base Position Qualities
  plot_qualities.py # Generates Quality_scores.png
  """
}

workflow {
    main:
        ch_reads = Channel.fromPath(params.reads_dir + "/*.fastq.gz") | collect
        
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

        MOTUS_PROFILE( BBDUK.out.trimmed_reads, params.meta_sample_name )

        FETCH_NCBI_GENOMES( 
          MOTUS_PROFILE.out.top_species | map { it.text }, 
          params.number_genome_representatives
          )
        
        BUILD_BOWTIE2_REFERENCE( 
          MOTUS_PROFILE.out.top_species | map { it.text }, 
          FETCH_NCBI_GENOMES.out.genomes
          )


        ALIGN_BOWTIE2( 
          BBDUK.out.trimmed_reads,
          BUILD_BOWTIE2_REFERENCE.out.reference,
          params.meta_sample_name
          )

        EXTRACT_PER_BASE_ALIGNMENT_SCORES(
          ALIGN_BOWTIE2.out.alignment,
          params.sample_count_alignments
          )

        EXTRACT_PER_BASE_PHRED(
          BBDUK.out.trimmed_reads
          )

        GENERATE_REPORTS(
          MEGAHIT.out.log, // For insert size mean and s.d.
          MULTIQC.out.data, // For number of molecules sequenced, for sequencing configuration, libary type, For plot, per position quality 
          RAW_FASTQC.out.fastqc_zip | map{ it[0]}, // For plot, per position quality
          EXTRACT_PER_BASE_ALIGNMENT_SCORES.out.table, // For plot, per postion alignment quality
          EXTRACT_PER_BASE_PHRED.out.table, // For plot, per postion alignment quality
          FETCH_NCBI_GENOMES.out.full_genomes_list, // For text report, species genome accessions
          MOTUS_PROFILE.out.top_species | map { it.text }, // For text report, top species
          ch_reads // For text report, index
          )
}