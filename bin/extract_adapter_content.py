#! /usr/bin/env python

import pandas as pd

df = pd.read_csv("adapter_content.txt", sep="\t", comment=">")

most_prevalent = df.describe().T['max'].idxmax()

# Categorize as Nextera or TruSeq
## Nextera is just `Nextera Transposase Sequence`
## TruSeq is anything the includes `Illumina`
## Other will be anything else

if "Nextera" in most_prevalent:
    adapter_type = "Nextera"
elif "Illumina" in most_prevalent:
    adapter_type = "TruSeq"
else:
    adapter_type = "Other"

print(f"Library Type (inferred by adapters in raw fastQC): {adapter_type}")