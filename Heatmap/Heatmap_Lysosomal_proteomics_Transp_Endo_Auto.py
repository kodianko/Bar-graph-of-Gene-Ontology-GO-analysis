import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# =========================
# Load data
# =========================
df1_path = r'X:\X\X\X\Lysosomal_proteomics_data.xlsx'
df1_tbl = pd.read_excel(df1_path)



plt.rcParams['font.family'] = 'Arial'

# =========================
# Groups
# =========================

groups = {
    'Vesicle_Mediated_Transport': ["RAB5C", "RBSN", "LRP1", "ADRB2", "AP2A2", "AP2B1", "AP2S1", "SORT1", "AP1S2", "AP1M2", "AP3D1", "AP3S1", "VPS54", "VPS4A", "VPS4B", "CHMP2A", "CHMP4A", "CHMP5", "CHMP6", "RAB7B", "NCOA4", "SQSTM1", "NEDD4", "RAB14", "RAB34", "HOOK1", "HOOK2", "SNX16", "CPNE1", "NPC1", "WDR81", "PLEKHM1", "DENND3", "RHOB", "TPCN2", "VPS41", "NSF", "VAMP7", "STX3", "DNAJC5", "DTX3L", "PLEKHF1", "ATP6V1H"],

    'Clathrin_Dependent_Endocytosis': ["CAV1", "NEDD4", "CLTC", "AP2S1", "CLTA", "AP2A1", "ADRB2", "AP2B1", "AP2A2", "LRP1"],

    'Autophagy': ["CLEC16A","FLCN","GNAI3","KLHL22","MLST8","NPRL2","RHEB","SH3BP4","SRC","TSC2","MAP1LC3B2","GABARAP","GABARAPL2","GABARAPL1","NEDD4","SQSTM1","CALCOCO2","NBR1","BNIP3","TECPR1","VPS36","CHMP2A","CHMP4A","CHMP5","CHMP6","VPS4A","VPS4B","WIPI1","WDR45","WDR45B","WDR81","FYCO1","UVRAG","ZFYVE26","SNX14","PLEKHM1","VPS41","EPG5","VAMP7","TPCN1","PLEKHF1","TMEM59","ATP6V1A","ATP6V1B2","ATP6V1H"],

}

# =========================
# Sample columns
# =========================

raw_sample_cols = [
    '7C1', '7C2', '7C3',
    '7N1', '7N2', '7N3',
    '7T1', '7T2', '7T3',
    '10C1', '10C2', '10C3',
    '10N1', '10N2', '10N3',
    '10T1', '10T2', '10T3'
]

sample_cols = [
    '7C1', '7C2', '7C3',
    '7T1', '7T2', '7T3',
    '10C1', '10C2', '10C3',
    '10T1', '10T2', '10T3',
    '7N1', '7N2', '7N3',
    '10N1', '10N2', '10N3'
]

# =========================
# Normalize columns
# =========================

for col in raw_sample_cols:
    df1_tbl[col + '_norm'] = df1_tbl[col] / df1_tbl[col].sum()

sample_cols_norm = [c + '_norm' for c in sample_cols]

# =========================
# Remove missing values
# =========================

non_nan = df1_tbl.dropna(subset=raw_sample_cols).copy()

# =========================
# Save folder
# =========================

save_to_folder = r'X:\X\X\X\X'
os.makedirs(save_to_folder, exist_ok=True)

norm_settings = {
    'Vesicle_Mediated_Transport': 1e-7,
    'Autophagy': 1e-7,
    'Clathrin_Dependent_Endocytosis': 1e-6,
    'Ferroptosis': 1e-6,
    'Exocytosis': 1e-6
}

# =========================
# Make heatmaps
# =========================

for group_name, gene_list in groups.items():

    df_subset = non_nan[
        (non_nan['Gene name'].isin(gene_list)) &
        (non_nan['Accession'].astype(str).str.strip() != 'P0DI83')
    ].copy()

    # Keep proteins upregulated in at least one comparison
    df_subset1 = df_subset[
        (
            (df_subset['7C/10C_pvalue'] < 0.05) &
            (df_subset['7C/10C_log2FC'] >= 1)
        ) |
        (
            (df_subset['7T/10T_pvalue'] < 0.05) &
            (df_subset['7T/10T_log2FC'] >= 1)
        ) |
        (
            (df_subset['7N/10N_pvalue'] < 0.05) &
            (df_subset['7N/10N_log2FC'] >= 1)
        )
    ].copy()

    df_subset_ordered = (
        df_subset1
        .set_index('Gene name')
        .reindex(gene_list)
        .dropna(how='all')
        .reset_index()
    )

    if df_subset_ordered.empty:
        print(f'No proteins found for {group_name}')
        continue

    matrix = df_subset_ordered[sample_cols_norm].copy()
    ytick = df_subset_ordered['Gene name']
    vmin = norm_settings[group_name]

    g = sns.clustermap(
        matrix,
        col_cluster=False,
        row_cluster=False,
        cmap='coolwarm',
        yticklabels=ytick,
        xticklabels=sample_cols,
        norm=LogNorm(
            vmin=vmin,
            vmax=matrix.max().max()
        ),
        figsize=(10, 12),
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

    g.fig.suptitle(group_name.replace('_', ' '), y=1.05)

    figurename = f'{group_name}_upregulated_heatmap'

    g.fig.savefig(
        os.path.join(save_to_folder, f'{figurename}.png'),
        dpi=300,
        bbox_inches='tight'
    )

    g.fig.savefig(
        os.path.join(save_to_folder, f'{figurename}.eps'),
        format='eps',
        bbox_inches='tight'
    )

    plt.show()
    plt.close(g.fig)

print('All heatmaps saved.')
