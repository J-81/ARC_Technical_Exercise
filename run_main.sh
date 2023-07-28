conda activate main

MAIN='/workspace/ARC_Technical_Exercise/main.nf'

mkdir /workspace/nftest
cd /workspace/nftest

nextflow run $MAIN \
    --reads_dir /workspace/data \
    --bbduk_adapters /workspace/bbmap_install/bbmap/resources/adapters.fa \
    --meta_sample_name ARC_TEST_SAMPLE \
    $@