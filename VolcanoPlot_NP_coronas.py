import bioinfokit
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import os
import seaborn as sns
from bioinfokit import analys, visuz

# =========================
# Load data
# =========================
df_path = r'X:\X\X\X\NP_coronas_data.xlsx'
df = pd.read_excel(df_path)


df = df.rename(columns={
    'Normalized1 FBS': 'FBS1',
    'Normalized2 FBS': 'FBS2',
    'Normalized3 FBS': 'FBS3',

    'Normalized1 TMA NPs': 'TMA1',
    'Normalized2 TMA NPs': 'TMA2',
    'Normalized3 TMA NPs': 'TMA3',

    'Normalized1 TMA:MUA MCNPs': 'MUA1',
    'Normalized2 TMA:MUA MCNPs': 'MUA2',
    'Normalized3 TMA:MUA MCNPs': 'MUA3',

    'Normalized1 TMA:MUS MCNPs': 'MUS1',
    'Normalized2 TMA:MUS MCNPs': 'MUS2',
    'Normalized3 TMA:MUS MCNPs': 'MUS3',

    'Normalized1 TMA:MUP MCNPs': 'MUP1',
    'Normalized2 TMA:MUP MCNPs': 'MUP2',
    'Normalized3 TMA:MUP MCNPs': 'MUP3'
})

# sample columns
sample_cols = [
    'FBS1', 'FBS2', 'FBS3',
    'TMA1', 'TMA2', 'TMA3',
    'MUA1', 'MUA2', 'MUA3',
    'MUS1', 'MUS2', 'MUS3',
    'MUP1', 'MUP2', 'MUP3'
]


# log2 transform
for col in sample_cols:
    df[f'log2({col})'] = np.log2(df[col] + 1e-6)

# log2 fold changes over FBS
# log2 fold changes over FBS
df['MUA/TMA_log2FC'] = (
    df[['log2(MUA1)', 'log2(MUA2)', 'log2(MUA3)']].mean(axis=1)
    -
    df[['log2(TMA1)', 'log2(TMA2)', 'log2(TMA3)']].mean(axis=1)
)

df['MUS/TMA_log2FC'] = (
    df[['log2(MUS1)', 'log2(MUS2)', 'log2(MUS3)']].mean(axis=1)
    -
    df[['log2(TMA1)', 'log2(TMA2)', 'log2(TMA3)']].mean(axis=1)
)

df['MUP/TMA_log2FC'] = (
    df[['log2(MUP1)', 'log2(MUP2)', 'log2(MUP3)']].mean(axis=1)
    -
    df[['log2(TMA1)', 'log2(TMA2)', 'log2(TMA3)']].mean(axis=1)
)


# WITH FDR
# Common proteins
common_MUA = df[
    (df['Relative mean TMA NPs'] > 0) &
    (df['Relative mean TMA:MUA MCNPs'] > 0)
]

common_MUS = df[
    (df['Relative mean TMA NPs'] > 0) &
    (df['Relative mean TMA:MUS MCNPs'] > 0)
]

common_MUP = df[
    (df['Relative mean TMA NPs'] > 0) &
    (df['Relative mean TMA:MUP MCNPs'] > 0)
]

plt.rcParams['font.family'] = 'Arial'

common_dfs = {
    'MUA': common_MUA,
    'MUS': common_MUS,
    'MUP': common_MUP
}

comparisons = {
    'MUA': ['MUA/TMA_log2FC', 'MUA vs TMA_Tukey_FDR'],
    'MUS': ['MUS/TMA_log2FC', 'MUS vs TMA_Tukey_FDR'],
    'MUP': ['MUP/TMA_log2FC', 'MUP vs TMA_Tukey_FDR']
}

for name in comparisons:

    fc_col, pval_col = comparisons[name]
    plot_df = common_dfs[name][[fc_col, pval_col, 'Gene Name']].dropna().copy()

    plot_df['highlight'] = (
        (plot_df[pval_col] < 0.05) &
        (abs(plot_df[fc_col]) >= 1)
    )

    highlighted = plot_df[plot_df['highlight']]

    plt.figure(figsize=(4.6, 6))

    plt.scatter(
        plot_df[fc_col],
        -np.log10(plot_df[pval_col]),
        s=8,
        color='lightgray'
    )

    plt.scatter(
        highlighted[fc_col],
        -np.log10(highlighted[pval_col]),
        s=12,
        color='#dd1c77',

    )

    plt.axvline(-1, color='grey', linestyle='--')
    plt.axvline(1, color='grey', linestyle='--')
    plt.axhline(-np.log10(0.05), color='grey', linestyle='--')

    plt.xlabel('log2FC')
    plt.ylabel('-log10(Adj. p-value)')
    plt.title(f'Common proteins: {name}/TMA')

    plt.xlim(-9, 12)
    plt.ylim(-0.5, 9)

    plt.legend()
    plt.tight_layout()
    # plt.show()

    save_folder = r'X:\X\X\X'

    plt.savefig(
        fr'{save_folder}\{name}_volcano_Adj.png',
        dpi=300,
        bbox_inches='tight'
    )

    plt.savefig(
        fr'{save_folder}\{name}_volcano_Adj.eps',
        format='eps',
        bbox_inches='tight'
    )

    plt.show()

print('gg')