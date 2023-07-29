eval "$(micromamba shell hook --shell bash)" \
    && micromamba activate main \
    && nextflow run main.nf \
        --bbduk_adapters /opt/conda/envs/main/bin/bbmap/resources/adapters.fa \
        --reads_dir test_data \
/workspace/ARC_Technical_Exercise/modules        --meta_sample_name ARC_TEST_SAMPLE \
        --number_genome_representatives 4 \
        --sample_count_alignments 10000 \
        --output final_test
