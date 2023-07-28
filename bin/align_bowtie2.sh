ORG_NAME="Bacteroides_eggerthii_1_reference"
ORG_NAME="Bacteroides_eggerthii"
ORG_NAME="Bacteroides_thetaiotaomicron"

bowtie2 -x ${ORG_NAME} \
    -1 '/workspace/data/read1.fastq.gz' \
    -2 '/workspace/data/read2.fastq.gz'  \
    --no-unal \
    -p 8 \
    -S SAMPLE_raw_${ORG_NAME}.sam

bowtie2 -x ${ORG_NAME} \
    -1 '/workspace/nftest/bbduk_out/read1-trimmed.fastq.gz' \
    -2 '/workspace/nftest/bbduk_out/read2-trimmed.fastq.gz'  \
    --no-unal \
    -p 8 \
    -S SAMPLE_trimmed_${ORG_NAME}.sam

