#! /usr/bin/env python

from Bio import SeqIO
import numpy as np
import pandas as pd
import gzip
import sys

INPUT_FILES=sys.argv[1]

# Initialize a list to hold quality scores for each base position
quality_scores_by_position = []

# Initialize maximum read length
max_len = 0

for read_file in INPUT_FILES.split(','):

    # TODO handle non-gzipped files
    # Open the gzipped FASTQ file
    with gzip.open(read_file, "rt") as handle:


        # Iterate over each record in the FASTQ file
        for record in SeqIO.parse(handle, "fastq"):

            # Convert Phred quality scores to integers and add them to our list
            quality_scores_by_position.append(record.letter_annotations["phred_quality"])

            # Update maximum read length
            max_len = max(max_len, len(record.letter_annotations["phred_quality"]))

# Pad shorter reads with np.nan to match longest read length
for i in range(len(quality_scores_by_position)):
    if len(quality_scores_by_position[i]) < max_len:
        quality_scores_by_position[i] += [np.nan] * (max_len - len(quality_scores_by_position[i]))

# Convert the list to a 2D numpy array
quality_scores_by_position = np.array(quality_scores_by_position)

# Compute the average quality scores for each base position
average_quality_scores_by_position = np.nanmean(quality_scores_by_position, axis=0)

# Create a DataFrame
df = pd.DataFrame({
    "Position": np.arange(1, len(average_quality_scores_by_position) + 1),
    "Average Quality Score": average_quality_scores_by_position,
})

# Print the DataFrame
df.to_csv("per_positon_quality_scores.csv", index=False)
