#! /usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

df_align = pd.read_csv("per_positon_alignment_scores.csv").set_index("Position")
df_qual = pd.read_csv("per_positon_quality_scores.csv").set_index("Position")

print(df_align.head())
print(df_qual.head())

df = df_align.join(df_qual)

print(df.head())

df['Quality Score'].plot(label='Alignment Inferred Quality Score')
df['Average Quality Score'].plot(label='Instrument Quality Score')

plt.ylabel('-10*log10(error rate)')
plt.legend(loc='best')

# Save the figure as a PNG file
plt.savefig('Quality_scores.png', dpi=600, bbox_inches='tight')