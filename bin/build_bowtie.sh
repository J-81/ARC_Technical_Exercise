ORG_NAME=$1

bowtie2-build genomes/*.fna.gz ${ORG_NAME// /_}
