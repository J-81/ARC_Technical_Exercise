bbduk.sh in=sample-1-R1-raw.fastq.gz in2=sample-1-R2-raw.fastq.gz out1=sample-1-R1-trimmed.fastq.gz \
         out2=sample-1-R2-trimmed.fastq.gz ref=ref-adapters.fa ktrim=l k=17 ftm=5 qtrim=rl \
         trimq=10 mlf=0.5 maxns=0 > bbduk.log 2>&1