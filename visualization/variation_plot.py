#!/usr/bin/env python3 
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df_gc = pd.read_csv('variation_endosymb_only_gc_genome.csv')
df_size = pd.read_csv('variation_endosymb_only_size.csv')
df_final = pd.DataFrame()
df_final['species'] = df_gc['species']
df_final['gc_std'] = df_gc['std']
df_final['size_std'] = df_size['std']

plt.figure(figsize=(8,8))
sns.scatterplot(data=df_final, x='gc_std', y='size_std')
plt.title('Standard Deviation of GC% and Size\nEndosymbionts Only', fontsize=14)
plt.xlabel('GC% Std. Deviation')
plt.ylabel('Size Std. Deviation')

# Only use this if size_std is actually in base pairs and large values
# plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f} Mb'))

plt.margins(x=0.1)

# Calculate quantiles correctly
gc_high, gc_low = df_final['gc_std'].quantile(0.65), df_final['gc_std'].quantile(0.4)
size_high, size_low = df_final['size_std'].quantile(0.65), df_final['size_std'].quantile(0.4)

select_high = df_final[(df_final['gc_std'] > gc_high) & (df_final['size_std'] > size_high)]
select_low = df_final[(df_final['gc_std'] < gc_low) & (df_final['size_std'] < size_low)]
size_extreme = df_final[df_final['size_std'] > df_final['size_std'].quantile(0.9)]
gc_extreme = df_final[df_final['gc_std'] > df_final['gc_std'].quantile(0.9)]

# Combine and remove duplicates for labeling
to_label = pd.concat([select_high, size_extreme, gc_extreme]).drop_duplicates()

for _, row in to_label.iterrows():
    plt.text(row['gc_std'] + 0.001, row['size_std'] + 0.001, row['species'], 
             fontsize=8, alpha=0.8)

print("Low variation species:")
print(select_low['species'].tolist())

plt.tight_layout()
plt.savefig('variation_plot.pdf', pad_inches=0.3, dpi=300)
plt.close()
