#! /usr/bin/env python

import sys

import pandas as pd

MULTIQC_DATA_TABLE=sys.argv[1]

df = pd.read_csv(MULTIQC_DATA_TABLE, sep="\t")

print(f"Sequenced Reads: {int(df['FastQC_mqc-generalstats-fastqc-total_sequences'].max())}")