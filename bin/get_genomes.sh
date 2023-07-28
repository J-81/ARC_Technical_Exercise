#! /usr/bin/env bash

ORG_NAME=$1 #  e.g. 'Bacteroides eggerthii'
LIMIT_TO=$2 #  e.g. 5

TAX_ID=$(esearch -db taxonomy -query "${ORG_NAME}" | \
    efetch -format xml | \
    xtract -pattern Taxon -element TaxId)

echo "Found taxonomy id ${TAX_ID} for ${ORG_NAME}"

ARRAY_ACCESSIONS=$(esearch -db taxonomy -query "txid${TAX_ID}[Organism:exp]" | elink -target assembly | \
    efetch -format docsum  | \
    xtract -pattern DocumentSummary -element RefSeq)

# Write the array to a file
printf "%s\n" ${ARRAY_ACCESSIONS[@]} > accession_ids.txt

# Initialize a counter
COUNTER=0

for A in $ARRAY_ACCESSIONS;
    do
        # If the counter reaches the limit, break the loop
        if [ $COUNTER -eq $LIMIT_TO ]; then
            break
        fi

        wget "$(esearch -db assembly -query ${A} | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}')"

        # Increment the counter
        ((COUNTER++))
    done