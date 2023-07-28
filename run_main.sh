conda activate main

MAIN='/workspace/ARC_Technical_Exercise/main.nf'

mkdir /workspace/nftest
cd /workspace/nftest

nextflow run $MAIN \
        -stub-run \
        --bbduk_adapters /workspace/bbmap_install/bbmap/resources/adapters.fa \
        --reads_dir /workspace/ARC_Technical_Exercise/test_data \
        --meta_sample_name ARC_TEST_SAMPLE \
        --stub_motus_results /workspace/ARC_Technical_Exercise/stub_data/taxonomy_profile.txt \
        --number_genome_representatives 4 \
        --sample_count_alignments 10000 \
        $@
