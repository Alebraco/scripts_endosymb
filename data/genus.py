import pandas as pd
import numpy as np

#Read bacterial metadata
df = pd.read_csv('genome_metadata.tsv', sep='\t', 
                 usecols=['Assembly Accession', 'Organism Name', 
                         'Assembly Stats Total Sequence Length',
                         'Annotation Count Gene Protein-coding',
                         'Annotation BUSCO Complete', 
                         'CheckM marker set', 'CheckM completeness',
                         'CheckM contamination'])
df = df.rename(columns={'Assembly Stats Total Sequence Length': 'Genome Size',})
# df = df[df['Assembly Accession'].str.contains('GCF_')]
# df = df[df['Organism Name'].str.contains('Rickettsia')]

#Retrieve endosymbionts if in the name
df['endosymb'] = df['Organism Name'].str.contains('endosymbiont')
endosymb_sp = df[df['endosymb']].copy()
endosymb_sp['Organism Name'] = endosymb_sp['Organism Name'].str.split('endosymbiont').str[0] + 'endosymbiont'

#Only retain those that have the following columns
df = df.dropna(subset=['Genome Size', 'Annotation Count Gene Protein-coding' 
                          ,'CheckM completeness', 'CheckM contamination'
                         ])
#Check if these genomes pass the CheckM filters
df = df[(df['CheckM completeness'] > 0.85) & (df['CheckM contamination'] < 0.10)]

#Clean the name of the organism
df['genus'] = (
    df['Organism Name']
    .str.replace(r"^(\'|\[|Candidatus |uncultured |unicellular |unidentified )*",'', regex=True) # prefixes
    .str.replace(r"[\]\']", '', regex=True) # suffixes
    .str.split().str[0]
)
#Define threshold for outliers (30% below the median value)
threshold = 0.3
crit = 1 # 1 for genome size, 0 for ratio

if crit == 1:
    # Genome size outliers
    df = df[['genus', 'Organism Name', 'Genome Size', 'endosymb']]
    df['median_size'] = df.groupby('genus')['Genome Size'].transform('median').round(0)
    df['outlier'] = (df['Genome Size'] < (df['median_size'] * (1 - threshold)))
    metrics = ['Genome Size', 'median_size']

else:
    # Ratio outliers
    df['gene_ratio'] = df['Annotation Count Gene Protein-coding'] / (df['Genome Size'] / 1000000)
    df['gene_ratio'] = df['gene_ratio'].round(2)

    df = df[['genus', 'Organism Name', 'gene_ratio', 'endosymb']]
    df['median_ratio'] = result.groupby('genus')['gene_ratio'].transform('median').round(2)
    df['outlier'] = df['gene_ratio'] < (df['median_ratio'] * (1 - threshold))
    metrics = ['gene_ratio', 'median_ratio']

endosymbs = endosymb_sp[['Organism Name', metrics[0]]].copy()
endosymbs[metrics[1]] = np.nan


outliers = df[df['outlier']][['Organism Name', metrics[0], metrics[1]]]
#Selects RECOGNIZED endosymbionts and OUTLIERS (based on criterion used)
symbs = pd.concat([endosymbs, outliers]).sort_values(by='Organism Name')

unique_symbs = symbs.drop_duplicates(subset=['Organism Name'])[['Organism Name']]
unique_symbs.to_csv('outliers.tsv', sep='\t', index=False)
