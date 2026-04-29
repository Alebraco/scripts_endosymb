#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
from delta_matrix import delta_matrix
from gcsize_dict import genome_gcsize
from sequence_features import collect_codon_stats
from utils import files_dir, genome_gcsize_json_path, load_or_compute, group_names

plot_dir = os.path.join('plots', 'igs_lengths')
os.makedirs(plot_dir, exist_ok=True)
group_colors = {
    'Endosymbionts' : '#FC8D62',
    'Free-Living Relatives' : '#66C2A5',
}

GROUP_RENAME = {
    'endosymb_only': 'Endosymbionts',
    'relatives_only': 'Free-Living Relatives',
}
GROUPS = list(group_colors.keys())


def stat_test(a, b):
    """Return (test_name, stat, pval) using Welch t-test if both groups are normal, else Mann-Whitney U."""
    p_sh_a = stats.shapiro(a).pvalue if 3 <= len(a) <= 5000 else None
    p_sh_b = stats.shapiro(b).pvalue if 3 <= len(b) <= 5000 else None
    if p_sh_a is not None and p_sh_b is not None and p_sh_a > 0.05 and p_sh_b > 0.05:
        stat, pval = stats.ttest_ind(a, b, equal_var=False, nan_policy='omit')
        return 'Welch t-test', stat, pval
    stat, pval = stats.mannwhitneyu(a, b, alternative='two-sided')
    return 'Mann-Whitney U', stat, pval


def plot_gene_counts():
    gene_counts_path = os.path.join(files_dir, 'gene_counts.csv')
    if not os.path.exists(gene_counts_path):
        print(f'Warning: gene_counts.csv not found at {gene_counts_path}. Skipping gene count boxplot.')
        return

    df = pd.read_csv(gene_counts_path)
    df['Group'] = df['Group'].replace(GROUP_RENAME)
    df = df[df['Group'].isin(GROUPS)]

    g_end = df.loc[df['Group'] == 'Endosymbionts', 'Gene_Count'].dropna()
    g_rel = df.loc[df['Group'] == 'Free-Living Relatives', 'Gene_Count'].dropna()
    test_name, stat, pval = stat_test(g_end, g_rel)
    print(f'{test_name} result: statistic={stat:.4f}, p-value={pval:.4g}')

    plt.figure(figsize=(6, 6))
    ax = sns.boxplot(data=df, x='Group', y='Gene_Count', order=GROUPS,
                     hue='Group', palette=group_colors, fliersize=0)
    ax.text(0.5, 0.95, f'{test_name}: p={pval:.4g}', transform=ax.transAxes,
            ha='center', va='top', fontsize=10)
    plt.xlabel('Group')
    plt.ylabel('Gene Count')
    plt.title('Gene Count by Group')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'gene_count_boxplot.pdf'))
    plt.close()


def plot_gene_lengths():
    gene_lengths_path = os.path.join(files_dir, 'mean_gene_lengths.csv')
    if not os.path.exists(gene_lengths_path):
        print(f'Warning: gene_lengths.csv not found at {gene_lengths_path}. Skipping gene length boxplot.')
        return

    df = pd.read_csv(gene_lengths_path)
    df['Group'] = df['Group'].replace(GROUP_RENAME)
    df = df[df['Group'].isin(GROUPS)]

    g_end = df.loc[df['Group'] == 'Endosymbionts', 'mean_gene_length'].dropna()
    g_rel = df.loc[df['Group'] == 'Free-Living Relatives', 'mean_gene_length'].dropna()
    test_name, stat, pval = stat_test(g_end, g_rel)
    print(f'{test_name} result: statistic={stat:.4f}, p-value={pval:.4g}')

    plt.figure(figsize=(6, 6))
    ax = sns.boxplot(data=df, x='Group', y='mean_gene_length', order=GROUPS,
                     hue='Group', palette=group_colors, fliersize=0)
    ax.text(0.5, 0.95, f'{test_name}: p={pval:.4g}', transform=ax.transAxes,
            ha='center', va='top', fontsize=10)
    plt.xlabel('Group')
    plt.ylabel('Gene Length')
    plt.title('Gene Length by Group')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'gene_length_boxplot.pdf'))
    plt.close()


def plot_igs_total_boxplot():
    df = pd.read_csv(os.path.join(files_dir, 'all_IGS_data.csv'))
    df['Group'] = df['Group'].replace(GROUP_RENAME)
    df = df[df['Group'].isin(GROUPS)]
    df_median = df.groupby(['Group', 'Species', 'File'])['IGS_Size'].mean().reset_index()

    plt.figure(figsize=(12, 12))
    sns.boxplot(data=df_median, x='IGS_Size', y='Species', hue='Group',
                palette=group_colors, hue_order=GROUPS, fliersize=2, dodge=True, gap=0.1)
    plt.xlabel('Species')
    plt.ylabel('IGS Size (bp)')
    plt.title('Intergenic Space Size Distribution by Species and Group')
    plt.legend(title='Group')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'IGS_total_boxplot.pdf'))
    plt.close()


def plot_igs_mean_boxplot(summary_df):
    plt.figure(figsize=(12, 12))
    sns.boxplot(data=summary_df, x='Group', y='mean_mean_IGS', hue='Group',
                palette=group_colors, hue_order=GROUPS, fliersize=0)
    sns.stripplot(data=summary_df, x='Group', y='mean_mean_IGS', color='black', alpha=0.7)
    plt.title('Intergenic Space (IGS) Size by Group', fontsize=24, fontweight='bold')
    plt.xticks(fontsize=18)
    plt.xlabel('Group', fontsize=20, fontweight='bold')
    plt.yticks(fontsize=18)
    plt.ylabel('IGS Size (bp)', fontsize=20, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'IGS_mean_boxplot.pdf'))
    plt.close()


def plot_igs_std_boxplot():
    genome_std_path = os.path.join(files_dir, 'genome_IGS.csv')
    if not os.path.exists(genome_std_path):
        print(f'Warning: genome_IGS.csv not found at {genome_std_path}. Skipping IGS std boxplot.')
        return

    df = pd.read_csv(genome_std_path)
    df['Group'] = df['Group'].replace(GROUP_RENAME)
    df = df[df['Group'].isin(GROUPS)]
    threshold = df['std_IGS'].quantile(0.90)
    df = df[df['std_IGS'] <= threshold]

    g_end = df.loc[df['Group'] == 'Endosymbionts', 'std_IGS'].dropna()
    g_rel = df.loc[df['Group'] == 'Free-Living Relatives', 'std_IGS'].dropna()
    test_name, stat, pval = stat_test(g_end, g_rel)
    print(f'IGS std {test_name} result: statistic={stat:.4f}, p-value={pval:.4g}')

    plt.figure(figsize=(12, 12))
    ax = sns.boxplot(data=df, x='Group', y='std_IGS', hue='Group',
                     palette=group_colors, hue_order=GROUPS, fliersize=0)
    # sns.stripplot(data=df, x='Group', y='std_IGS', hue='Group',
                #   palette=group_colors, hue_order=GROUPS, alpha=0.7, legend=False)
    ax.text(0.5, 0.95, f'{test_name}: p={pval:.4g}', transform=ax.transAxes,
            ha='center', va='top', fontsize=14)
    plt.title('Within-Genome IGS Standard Deviation by Group', fontsize=24, fontweight='bold')
    plt.xticks(fontsize=18)
    plt.xlabel('Group', fontsize=20, fontweight='bold')
    plt.yticks(fontsize=18)
    plt.ylabel('IGS Std Dev (bp)', fontsize=20, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'IGS_std_boxplot.pdf'))
    plt.close()


def plot_igs_std_species_boxplot():
    mean_path = os.path.join(files_dir, 'meanIGS.csv')
    if not os.path.exists(mean_path):
        print(f'Warning: meanIGS.csv not found at {mean_path}. Skipping IGS species std boxplot.')
        return

    df = pd.read_csv(mean_path)
    df['Group'] = df['Group'].replace(GROUP_RENAME)
    df = df[df['Group'].isin(GROUPS)]
    threshold = df['mean_std_IGS'].quantile(0.90)
    df = df[df['mean_std_IGS'] <= threshold]

    g_end = df.loc[df['Group'] == 'Endosymbionts', 'mean_std_IGS'].dropna()
    g_rel = df.loc[df['Group'] == 'Free-Living Relatives', 'mean_std_IGS'].dropna()
    test_name, stat, pval = stat_test(g_end, g_rel)
    print(f'IGS species std {test_name} result: statistic={stat:.4f}, p-value={pval:.4g}')

    plt.figure(figsize=(12, 12))
    ax = sns.boxplot(data=df, x='Group', y='mean_std_IGS', hue='Group',
                     palette=group_colors, hue_order=GROUPS, fliersize=0)
    ax.text(0.5, 0.95, f'{test_name}: p={pval:.4g}', transform=ax.transAxes,
            ha='center', va='top', fontsize=14)
    plt.title('Mean Within-Genome IGS Standard Deviation by Group\n(Per Species)', fontsize=24, fontweight='bold')
    plt.xticks(fontsize=18)
    plt.xlabel('Group', fontsize=20, fontweight='bold')
    plt.yticks(fontsize=18)
    plt.ylabel('Mean IGS Std Dev (bp)', fontsize=20, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'IGS_std_species_boxplot.pdf'))
    plt.close()


def plot_igs_scatterplot(summary_df):
    pivot_df = summary_df.pivot_table(index='Species', columns='Group', values='mean_mean_IGS').reset_index()

    plt.figure(figsize=(12, 12))
    sns.scatterplot(data=pivot_df, x='Endosymbionts', y='Free-Living Relatives',
                    s=100, color='gray', edgecolor='black')
    common_max = max(pivot_df['Endosymbionts'].max(), pivot_df['Free-Living Relatives'].max()) * 1.1
    plt.xlim(0, common_max)
    plt.ylim(0, common_max)
    plt.axline((0, 0), slope=1, color='gray', linestyle='--', linewidth=1, label='y=x (Equal IGS Size)')
    plt.xlabel('Endosymbionts', fontsize=20, fontweight='bold')
    plt.ylabel('Free-Living Relatives', fontsize=20, fontweight='bold')
    plt.margins(x=0.1)
    plt.title('Intergenic Space Size (bp)\nEndosymbionts vs Free-Living Relatives', fontsize=24, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'IGS_scatterplot.pdf'))
    plt.close()


def build_corr_df(summary_df):
    genome_json_path = genome_gcsize_json_path('endosymb+relatives')
    genome_data = load_or_compute(genome_json_path, genome_gcsize, 'endosymb+relatives')
    delta_df = delta_matrix(genome_data, mat_type='gc_genome')
    transposase = pd.read_csv(os.path.join(files_dir, 'median_transposases_species.csv'))
    gc4_df = pd.read_csv(os.path.join(files_dir, 'fourfold_gc_content.csv'))
    gc4_df = gc4_df[gc4_df['group'] == 'Endosymbionts']

    data_list = []
    for species, matrix in delta_df.items():
        sp_name = species.replace('_endosymbiont', '').replace('_', ' ')
        if sp_name not in summary_df['Species'].values:
            print(f'Skipping {species} as no IGS data is available.')
            continue
        print(f'Processing {species} for IGS vs Delta GC% correlation.')

        igs_size = summary_df[summary_df['Species'] == sp_name]['mean_mean_IGS'].values[0]

        ids = matrix.index.tolist()
        i, j = np.triu_indices_from(matrix, k=1)
        dist_list = []
        for index1, index2 in zip(i, j):
            id1, id2 = ids[index1], ids[index2]
            if ('_genomic' in id1) == ('_genomic' in id2):
                continue
            dist_list.append(matrix.loc[id1, id2])
        if not dist_list:
            print(f'No valid distance pairs for {species}, skipping.')
            continue

        if species not in genome_data:
            print(f'No genome GC% or size data for {species}, skipping.')
            continue

        gc_values = [genome_data[species].get(gid, {}).get('gc_genome')
                     for gid in genome_data[species]
                     if '_genomic' not in gid and isinstance(genome_data[species], dict)]
        size_values = [genome_data[species].get(gid, {}).get('size')
                       for gid in genome_data[species]
                       if '_genomic' not in gid and isinstance(genome_data[species], dict)]
        gc_value = np.median(gc_values) if gc_values else np.nan
        size = np.median(size_values) if size_values else np.nan

        if sp_name not in transposase['Species'].values:
            print(f'No transposase data for {sp_name}, skipping.')
            continue

        transposases = transposase[transposase['Species'] == sp_name]['Median_Transposases'].values[0]
        delta_gc_value = np.median(dist_list)
        gc4_row = gc4_df[gc4_df['species'] == sp_name]
        gc4_value = gc4_row['overall_gc4'].iloc[0] if not gc4_row.empty else np.nan
        data_list.append({'species': sp_name, 'mean_IGS': igs_size, 'DeltaGC': delta_gc_value,
                          'genome_gc': gc_value, 'genome_size': size, 'transposases': transposases,
                          'gc4': gc4_value})

    return pd.DataFrame(data_list)


CORR_X_COLS = ['DeltaGC', 'gc4', 'genome_gc', 'genome_size', 'transposases']

CORR_LABELS = {
    'species': {
        'DeltaGC':     'Delta GC% (median |Δ| over pairs)',
        'gc4':         'GC4 (4-fold-degenerate site GC%)',
        'genome_gc':   'Genome GC%',
        'genome_size': 'Genome Size (bp)',
        'transposases':'Median Number of Transposases',
    },
    'genome': {
        'DeltaGC':     'Delta GC% (|genome GC − median relative GC|)',
        'gc4':         'GC4 (4-fold-degenerate site GC%)',
        'genome_gc':   'Genome GC%',
        'genome_size': 'Genome Size (bp)',
        'transposases':'Total Transposases',
    },
}


def plot_igs_correlation_panel(corr_df, level, output_name):
    """
    Render all IGS-vs-metric scatterplots in a single 2x3 panel PDF.
    `level` is 'species' or 'genome' and selects axis labels / panel title.
    """
    labels = CORR_LABELS[level]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes_flat = axes.flatten()

    for ax, x_col in zip(axes_flat, CORR_X_COLS):
        sub = corr_df.dropna(subset=[x_col, 'mean_IGS'])
        sns.scatterplot(data=sub, x=x_col, y='mean_IGS', alpha=0.6, ax=ax)
        ax.set_xlabel(labels[x_col])
        ax.set_ylabel('Mean IGS Size (bp)')
        if len(sub) >= 3:
            rho, pval = stats.spearmanr(sub[x_col], sub['mean_IGS'])
            stat_str = f'Spearman ρ = {rho:.2f}, p = {pval:.2g} (n={len(sub)})'
        else:
            stat_str = f'Spearman ρ = n/a (n={len(sub)})'
        ax.set_title(f'{labels[x_col]}\n{stat_str}', fontsize=10)

    # Hide any unused axes (we have 5 metrics in a 6-cell grid)
    for ax in axes_flat[len(CORR_X_COLS):]:
        ax.set_visible(False)

    suptitle_suffix = 'per species' if level == 'species' else 'per genome'
    fig.suptitle(f'IGS Size vs Genomic Metrics in Endosymbionts ({suptitle_suffix})',
                 fontsize=16, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(plot_dir, output_name))
    plt.close(fig)


def build_genome_corr_df():
    """
    Build a per-endosymbiont-genome dataframe with one row per genome.
    Columns: species, file, mean_IGS, gc4, DeltaGC, genome_gc, genome_size, transposases.
    DeltaGC = |genome_gc - median(GC of matched relatives in same species group)|.
    GC4 is computed per genome from analysis/sequence_features.collect_codon_stats.
    Transposases come from files/transposase_summary.csv (Total_Transposases per File).
    """
    # Per-genome IGS
    igs_df = pd.read_csv(os.path.join(files_dir, 'genome_IGS.csv'))
    igs_df = igs_df[igs_df['Group'] == 'endosymb_only'].copy()

    # Per-accession genome GC, size + matched-relatives median GC
    genome_data = load_or_compute(
        genome_gcsize_json_path('endosymb+relatives'),
        genome_gcsize, 'endosymb+relatives'
    )
    rel_median_by_species = {}      # sp_name -> median GC of relatives
    endo_gc_by_file = {}            # file -> endosymb genome GC
    endo_size_by_file = {}          # file -> endosymb genome size
    for species, accessions in genome_data.items():
        sp_name = species.replace('_endosymbiont', '').replace('_', ' ')
        if not isinstance(accessions, dict):
            continue
        rel_gcs = [a.get('gc_genome') for gid, a in accessions.items()
                   if '_genomic' in gid and isinstance(a, dict)]
        rel_gcs = [g for g in rel_gcs if g is not None]
        if rel_gcs:
            rel_median_by_species[sp_name] = float(np.median(rel_gcs))
        for gid, a in accessions.items():
            if '_genomic' not in gid and isinstance(a, dict):
                if a.get('gc_genome') is not None:
                    endo_gc_by_file[gid] = a['gc_genome']
                if a.get('size') is not None:
                    endo_size_by_file[gid] = a['size']

    # Per-genome GC4 (cached because it walks every GFF + FASTA)
    gc4_cache = os.path.join(files_dir, 'codon_stats_endosymb_only.csv')
    if os.path.isfile(gc4_cache):
        codon_df = pd.read_csv(gc4_cache)
    else:
        codon_df = collect_codon_stats('endosymb_only', group_label='endosymb_only')
        os.makedirs(files_dir, exist_ok=True)
        codon_df.to_csv(gc4_cache, index=False)
    gc4_by_file = dict(zip(codon_df['File'], codon_df['GC4']))

    # Per-genome transposase counts
    transposase_by_file = {}
    transposase_path = os.path.join(files_dir, 'transposase_summary.csv')
    if os.path.exists(transposase_path):
        tdf = pd.read_csv(transposase_path)
        tdf = tdf[tdf['Group'] == 'endosymb_only']
        transposase_by_file = dict(zip(tdf['File'], tdf['Total_Transposases']))
    else:
        print(f'Warning: {transposase_path} not found; transposase column will be NaN.')

    rows = []
    for _, row in igs_df.iterrows():
        sp_name = row['Species']
        f = row['File']
        endo_gc = endo_gc_by_file.get(f)
        rel_med = rel_median_by_species.get(sp_name)
        delta_gc = abs(endo_gc - rel_med) if (endo_gc is not None and rel_med is not None) else np.nan
        rows.append({
            'species': sp_name,
            'file': f,
            'mean_IGS': row['mean_IGS'],
            'gc4': gc4_by_file.get(f, np.nan),
            'DeltaGC': delta_gc,
            'genome_gc': endo_gc if endo_gc is not None else np.nan,
            'genome_size': endo_size_by_file.get(f, np.nan),
            'transposases': transposase_by_file.get(f, np.nan),
        })
    return pd.DataFrame(rows)


def load_summary_df():
    df = pd.read_csv(os.path.join(files_dir, 'meanIGS.csv'))
    df['Group'] = df['Group'].replace(GROUP_RENAME)
    df = df[df['Group'].isin(GROUPS)]
    threshold = df['mean_mean_IGS'].quantile(0.90)
    return df[df['mean_mean_IGS'] <= threshold]


if __name__ == '__main__':
    plot_gene_counts()
    plot_gene_lengths()
    plot_igs_total_boxplot()
    summary_df = load_summary_df()
    plot_igs_mean_boxplot(summary_df)
    plot_igs_std_boxplot()
    plot_igs_std_species_boxplot()
    plot_igs_scatterplot(summary_df)
    corr_df = build_corr_df(summary_df)
    plot_igs_correlation_panel(corr_df, level='species',
                               output_name='IGS_correlations_endosymb_species.pdf')
    genome_corr_df = build_genome_corr_df()
    plot_igs_correlation_panel(genome_corr_df, level='genome',
                               output_name='IGS_correlations_endosymb_genome.pdf')
