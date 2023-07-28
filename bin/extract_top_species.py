#! /usr/bin/env python

import pandas as pd

df = pd.read_csv("taxonomy_profile.txt", 
                  comment="#", 
                  sep="\t",
                  names=["Species","Relative Abundance"]).sort_values(by="Relative Abundance", ascending=False)

# Extract and print top species by name
top_species = df["Species"].iloc[0].split('[')[0].strip()

with open("top_species.txt", "w") as f:
    f.write(top_species)