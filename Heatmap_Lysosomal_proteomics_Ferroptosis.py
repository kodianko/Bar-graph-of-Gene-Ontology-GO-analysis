import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =========================
# Load data
# =========================
df1_path = r'X:\X\X\Lysosomal_proteomics_data.xlsx'
df1_tbl = pd.read_excel(df1_path)

# group1 = ["SNX16", "SORT1", "NEDD4", "NCOA4", "PLEKHM1", "RBSN", "CHMP5", "NSF", "RAB5C", "VPS4A", "AP3D1", "AP2A1", "ADRB2", "AP2B1", 'AP2A2', "VAMP7", "RAB14", "AP1S2", "CPNE1", "AP2S1", "EPG5", "AP3S1", "STX3", "SQSTM1", "NPC1", "LRP1", "DENND3", "RHOB", "DNAJC5", "RAB34", "VPS41", "VPS53", "VPS54", "CHMP6", "TPCN2", "DTX3L", "HOOK2", "HOOK1", "VPS4B", "CHMP2A", "ATP6V1H", "AP1M2"]
# group1 = ["RAB5C", "RBSN", "LRP1", "ADRB2", "AP2A2", "AP2B1", "AP2S1", "SORT1", "AP1S2", "AP1M2", "AP3D1", "VPS54", "VPS4A", "VPS4B", "CHMP2A", "CHMP4A", "CHMP5", "CHMP6", "RAB7B", "NCOA4", "SQSTM1", "NEDD4", "RAB14", "RAB34", "HOOK1", "HOOK2", "SNX16", "CPNE1", "NPC1", "WDR81", "PLEKHM1", "DENND3", "RHOB", "TPCN2", "ANXA2", "VPS41", "NSF", "VAMP7", "STX3", "DNAJC5", "DTX3L", "PLEKHF1", "ATP6V1H"]
group1 = ["NSF", "SYNGR1", "VAMP7", "STX3", "VPS4A", "UNC13D", "DNAJC5", "SDC4", "RAB3D", "VPS4B", "CHMP2A", "SDC1"]


colnorm = df1_tbl[['7C1', '7C2', '7C3', '7N1', '7N2', '7N3', '7T1', '7T2', '7T3', '10C1', '10C2', '10C3', '10N1', '10N2', '10N3', '10T1', '10T2', '10T3']]
for _col in colnorm:
    df1_tbl[_col + '_norm'] = df1_tbl[_col]/df1_tbl[_col].sum()

non_nan = df1_tbl.dropna(subset=['7C1','7C2', '7C3', '7N1', '7N2', '7N3', '7T1', '7T2', '7T3', '10C1', '10C2', '10C3', '10N1', '10N2', '10N3', '10T1', '10T2', '10T3'])

df_subset = non_nan[
    (non_nan['Gene name'].isin(group1)) &
    (non_nan['Accession'].astype(str).str.strip() != 'P0DI83')
]

df_subset1 = df_subset[
    (
        (df_subset['7C/10C_pvalue'] < 0.05) & (abs(df_subset['7C/10C_log2FC']) >= 1)
    ) |
    (
        (df_subset['7T/10T_pvalue'] < 0.05) & (abs(df_subset['7T/10T_log2FC']) >= 1)
    ) |
    (
        (df_subset['7N/10N_pvalue'] < 0.05) & (abs(df_subset['7N/10N_log2FC']) >= 1)
    )
]

df_subset_ordered = (
    df_subset1
    .set_index('Gene name')
    .reindex(group1)
    .dropna(how='all')   # removes genes not present
    .reset_index()
)

sample_cols = [
    '7C1', '7C2', '7C3',
    '7T1', '7T2', '7T3',
    '10C1', '10C2', '10C3',
    '10T1', '10T2', '10T3',
    '7N1', '7N2', '7N3',
    '10N1', '10N2', '10N3'
]

sample_cols_norm = [c + '_norm' for c in sample_cols]

matrix = df_subset_ordered[sample_cols_norm].copy()
ytick = df_subset_ordered['Gene name']


g = sns.clustermap(
    matrix,
    col_cluster=False,
    row_cluster=False,
    cmap='coolwarm',
    yticklabels=ytick,
    xticklabels=sample_cols,
    norm=LogNorm(vmin=matrix.min().min(), vmax=matrix.max().max()),
    figsize=(10, 12),
    cbar_pos=[.28, .98, .625, .02],
    cbar_kws={'orientation': "horizontal", 'label': 'Normalized Abundance ($log_{10}$)'}
)

hm = g.ax_heatmap.get_position()
g.ax_heatmap.set_position([hm.x0*1.33, hm.y0*0.90, hm.width*0.9, hm.height*1.18])

# plt.show()



save_to_folder = r'E:\!!!Lysosomal Proteomics\Lysosomal proteins\figures\Heatmap\grouped'
figurename = "Exo_MCF7upandMCF10Adown_2026004228"
for extension in ['png', 'eps']:
    g.fig.savefig(f'{save_to_folder}/{figurename}.{extension}', dpi=300)


#FOR MCF10A
df_subset2 = df_subset[
    (
        (df_subset['7C/10C_pvalue'] < 0.05) & (df_subset['7C/10C_log2FC'] <= -1)
    ) |
    (
        (df_subset['7T/10T_pvalue'] < 0.05) & (df_subset['7T/10T_log2FC'] <= -1)
    ) |
    (
        (df_subset['7N/10N_pvalue'] < 0.05) & (df_subset['7N/10N_log2FC'] <= -1)
    )
]

df_subset_ordered2 = (
        df_subset2
        .set_index('Gene name')
        .reindex(group1)
        .dropna(how='all')  # removes genes not present
        .reset_index()
)

sample_cols2 = [
        '7C1', '7C2', '7C3',
        '7T1', '7T2', '7T3',
        '10C1', '10C2', '10C3',
        '10T1', '10T2', '10T3',
        '7N1', '7N2', '7N3',
        '10N1', '10N2', '10N3'
]

sample_cols_norm2 = [c + '_norm' for c in sample_cols2]

matrix2 = df_subset_ordered2[sample_cols_norm2].copy()
ytick2 = df_subset_ordered2['Gene name']

g = sns.clustermap(
        matrix2,
        col_cluster=False,
        row_cluster=False,
        cmap='coolwarm',
        yticklabels=ytick2,
        xticklabels=sample_cols2,
        norm=LogNorm(vmin=matrix2.min().min(), vmax=matrix2.max().max()),
        figsize=(10, 12),
        cbar_pos=[.28, .98, .625, .02],
        cbar_kws={'orientation': "horizontal", 'label': 'Normalized Abundance ($log_{10}$)'}
)

hm = g.ax_heatmap.get_position()
g.ax_heatmap.set_position([hm.x0 * 1.33, hm.y0 * 0.90, hm.width * 0.9, hm.height * 1.18])

plt.show()

save_to_folder = r'X:\X\X\X\X'
figurename = "Exo_MCF10Adown_2026004228"
for extension in ['png', 'eps']:
    g.fig.savefig(f'{save_to_folder}/{figurename}.{extension}', dpi=300)