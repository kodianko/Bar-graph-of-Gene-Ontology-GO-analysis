import os
import numpy as np
import pandas as pd
import matplotlib.pylab as plt

# =========================
# Load data
# =========================

go_path = r'X:\X\X\Lysosomal_proteomics_data.xlsx'
go_tbl = pd.read_excel(go_path)

sample_cols = [
    '7C1', '7C2', '7C3',
    '7N1', '7N2', '7N3',
    '7T1', '7T2', '7T3',
    '10C1', '10C2', '10C3',
    '10N1', '10N2', '10N3',
    '10T1', '10T2', '10T3'
]

go_tbl_nonan = go_tbl.dropna(subset=sample_cols)

df_all = go_tbl_nonan.loc[
    go_tbl_nonan['Accession'].astype(str).str.strip() != 'P0DI83'
].copy()

df_lyso = go_tbl_nonan.loc[
    (go_tbl_nonan['KnownLysosomalProteins'] == True) &
    (go_tbl_nonan['Accession'].astype(str).str.strip() != 'P0DI83')
].copy()

# =========================
# Groups
# =========================

groups = {
    'Ferroptosis': ['CYB561A3', 'LTF', 'MCOLN1', 'SLC11A2', 'SLC40A1', 'STEAP3', 'TFRC', 'FTH1', 'FTL', 'NCOA4', 'RAB7A', 'CTH', 'G6PD', 'GCLC', 'GCLM', 'GPX4', 'GSR', 'GSS', 'SLC3A2', 'SLC7A11', 'GLS', 'GLUD1', 'GOT1', 'GOT2', 'SLC1A5', 'SLC38A1', 'ACSL4', 'ALOX12', 'ALOX15', 'ALOX15B', 'LPCAT3', 'PEBP1', 'POR', 'ACACA', 'ACLY', 'ACSF2', 'FASN', 'FDFT1', 'SQLE', 'AIFM2', 'DHFR', 'DHODH', 'GCH1', 'AKR1C2', 'AKR1C3', 'PLA2G6', 'PRDX1', 'PRDX5', 'PRDX6', 'TXN', 'HMOX1', 'KEAP1', 'NFE2L2', 'NQO1', 'SRXN1', 'TXNRD1'],

    'Clathrin_Dependent_Endocytosis': ['CAV1', 'NEDD4', 'CLTC', 'AP2S1', 'CLTA', 'AP2A1', 'ADRB2', 'AP2B1', 'AP2A2', 'LRP1'],

    'Vesicle_Mediated_Transport': ['RAB5C', 'RBSN', 'LRP1', 'ADRB2', 'AP2A2', 'AP2B1', 'AP2S1', 'SORT1', 'AP1S2', 'AP1M2', 'AP3D1', 'AP3S1', 'VPS54', 'VPS4A', 'VPS4B', 'CHMP2A', 'CHMP4A', 'CHMP5', 'CHMP6', 'RAB7B', 'NCOA4', 'SQSTM1', 'NEDD4', 'RAB14', 'RAB34', 'HOOK1', 'HOOK2', 'SNX16', 'CPNE1', 'NPC1', 'WDR81', 'PLEKHM1', 'DENND3', 'RHOB', 'TPCN2', 'VPS41', 'NSF', 'VAMP7', 'STX3', 'DNAJC5', 'DTX3L', 'PLEKHF1', 'ATP6V1H'],

    'Autophagy': ['CLEC16A', 'FLCN', 'GNAI3', 'KLHL22', 'MLST8', 'NPRL2', 'RHEB', 'SH3BP4', 'SRC', 'TSC2', 'MAP1LC3B2', 'GABARAP', 'GABARAPL2', 'GABARAPL1', 'NEDD4', 'SQSTM1', 'CALCOCO2', 'NBR1', 'BNIP3', 'TECPR1', 'VPS36', 'CHMP2A', 'CHMP4A', 'CHMP5', 'CHMP6', 'VPS4A', 'VPS4B', 'WIPI1', 'WDR45', 'WDR45B', 'WDR81', 'FYCO1', 'UVRAG', 'ZFYVE26', 'SNX14', 'PLEKHM1', 'VPS41', 'EPG5', 'VAMP7', 'TPCN1', 'PLEKHF1', 'TMEM59', 'ATP6V1A', 'ATP6V1B2', 'ATP6V1H'],

    'Exocytosis': ['NSF', 'SYNGR1', 'VAMP7', 'STX3', 'VPS4A', 'UNC13D', 'DNAJC5', 'SDC4', 'RAB3D', 'VPS4B', 'CHMP2A', 'SDC1']
}

# =========================
# Comparisons
# =========================

comparisons = {
    '7N10N': {
        'title': '7N/10N',
        'fc_col': '7N/10N_log2FC',
        'p_col': '7N/10N_pvalue'
    },
    '7T10T': {
        'title': '7T/10T',
        'fc_col': '7T/10T_log2FC',
        'p_col': '7T/10T_pvalue'
    },
    '7C10C': {
        'title': '7C/10C',
        'fc_col': '7C/10C_log2FC',
        'p_col': '7C/10C_pvalue'
    }
}

# =========================
# Output folder
# =========================

save_to_folder = r'X:\X\X\X'
os.makedirs(save_to_folder, exist_ok=True)

plt.rcParams['font.family'] = 'Arial'

# =========================
# Volcano plots
# =========================

for comparison_name, comparison_info in comparisons.items():

    fc_col = comparison_info['fc_col']
    p_col = comparison_info['p_col']
    title = comparison_info['title']

    for group_name, gene_list in groups.items():

        if group_name == 'Ferroptosis':
            plot_df = df_all[
                df_all['Gene name'].isin(gene_list)
            ].dropna(subset=[fc_col, p_col]).copy()

            upregulated = plot_df[
                (plot_df[fc_col] >= 1) &
                (plot_df[p_col] < 0.05)
                ].copy()

            downregulated = plot_df[
                (plot_df[fc_col] <= -1) &
                (plot_df[p_col] < 0.05)
                ].copy()

            highlight_label = 'Upregulated ferroptosis-associated proteins'
            plot_title = f'{title}: Ferroptosis'

        else:
            plot_df = df_lyso.dropna(subset=[fc_col, p_col]).copy()

            highlighted = plot_df[
                plot_df['Gene name'].isin(gene_list) &
                (plot_df[p_col] < 0.05)
            ].copy()

            highlight_label = group_name.replace('_', ' ')
            plot_title = f'{title}: {group_name.replace("_", " ")}'

        plt.figure(figsize=(4.6, 6))

        # Background
        # For Ferroptosis: grey = 56 ferroptosis-associated proteins
        # For other groups: grey = known lysosomal proteins
        plt.scatter(
            plot_df[fc_col],
            -np.log10(plot_df[p_col]),
            s=5,
            color='lightgray',
            label='Ferroptosis-associated proteins' if group_name == 'Ferroptosis' else 'Known lysosomal proteins'
        )

        # Highlighted proteins
        plt.scatter(
            upregulated[fc_col],
            -np.log10(upregulated[p_col]),
            s=10,
            color='#dd1c77',
            label='Upregulated ferroptosis-associated proteins'
        )

        plt.scatter(
            downregulated[fc_col],
            -np.log10(downregulated[p_col]),
            s=10,
            color='#dd1c77',
            label='Downregulated ferroptosis-associated proteins'
        )

        plt.axvline(-1, color='grey', linestyle='--')
        plt.axvline(1, color='grey', linestyle='--')
        plt.axhline(-np.log10(0.05), color='grey', linestyle='--')

        plt.xlim(-9, 12)
        plt.ylim(-0.5, 10)

        plt.xlabel('log2FC')
        plt.ylabel('-log10(p-value)')
        plt.title(plot_title)
        plt.legend(fontsize=7)

        for _, row in upregulated.iterrows():
            plt.text(
                row[fc_col],
                -np.log10(row[p_col]),
                row['Gene name'],
                fontsize=7
            )

        for _, row in downregulated.iterrows():
            plt.text(
                row[fc_col],
                -np.log10(row[p_col]),
                row['Gene name'],
                fontsize=7
            )

        plt.tight_layout()

        figurename = f'{comparison_name}_{group_name}_volcano'

        plt.savefig(
            os.path.join(save_to_folder, f'{figurename}.png'),
            dpi=300,
            bbox_inches='tight'
        )

        plt.savefig(
            os.path.join(save_to_folder, f'{figurename}.eps'),
            format='eps',
            bbox_inches='tight'
        )

        plt.show()
        plt.close()

print('All volcano plots saved.')