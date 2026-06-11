import bioinfokit
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
from bioinfokit import analys, visuz
from matplotlib.colors import LogNorm
import os

# =========================
# Load main file only
# =========================
go_path = r'X:\X\X\Lysosomal_proteomics_data.xlsx'
df1_tbl = pd.read_excel(go_path)

# =========================
# Exocytosis / secretion protein list
# =========================
group1 = [
    "NSF", "SYNGR1", "VAMP7", "STX3", "VPS4A", "UNC13D",
    "DNAJC5", "SDC4", "RAB3D", "VPS4B", "CHMP2A", "SDC1"
]

# =========================
# Clean gene names
# =========================
df1_tbl['Gene_clean'] = (
    df1_tbl['Gene name']
    .astype(str)
    .str.upper()
    .str.strip()
)

df1_tbl['Group1'] = df1_tbl['Gene_clean'].isin(group1)

# =========================
# Sample columns
# =========================
sample_cols = [
    '7C1', '7C2', '7C3',
    '7T1', '7T2', '7T3',
    '10C1', '10C2', '10C3',
    '10T1', '10T2', '10T3',
    '7N1', '7N2', '7N3',
    '10N1', '10N2', '10N3'
]

# =========================
# Normalize abundance columns
# =========================
for col in sample_cols:
    df1_tbl[col + '_norm'] = df1_tbl[col] / df1_tbl[col].sum()

sample_cols_norm = [c + '_norm' for c in sample_cols]

# =========================
# Remove rows with missing values
# =========================
non_nan = df1_tbl.dropna(subset=sample_cols)

df = non_nan.loc[
    non_nan['Group1'] == True
].copy()

# =========================
# MCF7 significant changed proteins
# abs(log2FC) >= 1
# =========================
df_subset_mcf7 = df[
    (
        (df['7C/10C_pvalue'] < 0.05) &
        (abs(df['7C/10C_log2FC']) >= 1)
    ) |
    (
        (df['7T/10T_pvalue'] < 0.05) &
        (abs(df['7T/10T_log2FC']) >= 1)
    ) |
    (
        (df['7N/10N_pvalue'] < 0.05) &
        (abs(df['7N/10N_log2FC']) >= 1)
    )
].copy()

df_subset_mcf7['Heatmap_group'] = 'MCF7 up significant'

# =========================
# MCF10A downregulated proteins
# log2FC <= -1
# =========================
df_subset_mcf10a = df[
    (
        (df['7C/10C_pvalue'] < 0.05) &
        (df['7C/10C_log2FC'] <= -1)
    ) |
    (
        (df['7T/10T_pvalue'] < 0.05) &
        (df['7T/10T_log2FC'] <= -1)
    ) |
    (
        (df['7N/10N_pvalue'] < 0.05) &
        (df['7N/10N_log2FC'] <= -1)
    )
].copy()

df_subset_mcf10a['Heatmap_group'] = 'MCF10A down significant'

# =========================
# Combine both groups into one heatmap
# =========================
df_combined = pd.concat(
    [df_subset_mcf7, df_subset_mcf10a],
    axis=0
)

df_combined = df_combined.drop_duplicates(subset='Gene name')

# =========================
# Order proteins by custom list
# =========================
protein_order = group1

order_dict = {
    gene.upper(): i
    for i, gene in enumerate(protein_order)
}

df_combined['Gene_clean'] = (
    df_combined['Gene name']
    .astype(str)
    .str.upper()
    .str.strip()
)

df_combined['protein_order'] = df_combined['Gene_clean'].map(order_dict)

df_combined = df_combined[
    df_combined['Gene_clean'].isin([g.upper() for g in protein_order])
]

group_order = [
    'MCF7 significant',
    'MCF10A down'
]

df_combined['Heatmap_group'] = pd.Categorical(
    df_combined['Heatmap_group'],
    categories=group_order,
    ordered=True
)

df_combined = df_combined.sort_values(
    by=['Heatmap_group', 'protein_order']
)

# =========================
# Matrix for heatmap
# =========================
matrix = df_combined[sample_cols_norm].copy()
matrix = matrix.replace(0, np.nan)

ytick = df_combined['Gene name']

# =========================
# Heatmap
# =========================
g = sns.clustermap(
    matrix,
    col_cluster=False,
    row_cluster=False,
    cmap='coolwarm',
    yticklabels=ytick,
    xticklabels=sample_cols,
    norm=LogNorm(
        vmin=matrix.min().min(),
        vmax=matrix.max().max()
    ),
    figsize=(10, 8),
    cbar_pos=[.28, .98, .625, .02],
    cbar_kws={
        'orientation': "horizontal",
        'label': 'Normalized Abundance ($log_{10}$)'
    }
)

hm = g.ax_heatmap.get_position()
g.ax_heatmap.set_position([
    hm.x0 * 1.33,
    hm.y0 * 0.90,
    hm.width * 0.9,
    hm.height * 1.18
])

plt.show()

# =========================
# Save figure
# =========================
save_to_folder = r'X:\X\X\X\X\X'
os.makedirs(save_to_folder, exist_ok=True)

figurename = "Exo_MCF7significant_MCF10Adown_combined_20260611"

for extension in ['png', 'eps']:
    g.fig.savefig(
        f'{save_to_folder}/{figurename}.{extension}',
        dpi=300,
        bbox_inches='tight'
    )

print('Saved combined heatmap')
print('Number of proteins in heatmap:', len(df_combined))
print(df_combined[['Gene name', 'Heatmap_group']])