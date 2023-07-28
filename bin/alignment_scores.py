#! /usr/bin/env python

import pandas as pd
import pysam
import sys

SAM_INPUT=sys.argv[1]
SUBSAMPLE=int(sys.argv[2])

# Open the SAM file
file = pysam.AlignmentFile(SAM_INPUT,"r")

# Create an empty list to hold our data
data = []

# Iterate over each read in the file
for i, read in enumerate(file):
    if i % 1000 == 0 and i > 0:
        if SUBSAMPLE > 0:
            print(f"Processing alignment {i:6} of {SUBSAMPLE} target sample size.", end="\r")
        else:
            print(f"Processing alignment {i:6}.", end="\r")
    if i > SUBSAMPLE and SUBSAMPLE > 0:
        break
    # Get the query name, sequence and quality scores of the read
    query_name = read.query_name
    sequence = read.query_sequence
    quality_scores = read.query_qualities

    # Append this data to our list
    for position, (base, score) in enumerate(zip(sequence, quality_scores)):
        # Postion is 0-based, so we add 1 to it
        data.append([query_name, position+1, base, score])

# Convert our data to a DataFrame
df = pd.DataFrame(data, columns=['Read', 'Position', 'Base', 'Quality Score'])

# Close the file
file.close()

df.groupby("Position").mean('Quality Score').to_csv("per_positon_alignment_scores.csv")
