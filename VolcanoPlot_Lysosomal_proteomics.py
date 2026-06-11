import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import os

go_path = r'XE:\X\X\X\Lysosomal_proteomics_data.xlsx'
go_tbl = pd.read_excel(go_path)

# remove rows with missing values in sample columns
sample_cols = [
    '7C1', '7C2', '7C3',
    '7N1', '7N2', '7N3',
    '7T1', '7T2', '7T3',
    '10C1', '10C2', '10C3',
    '10N1', '10N2', '10N3',
    '10T1', '10T2', '10T3'
]

go_tbl_nonan = go_tbl.dropna(subset=sample_cols)

# keep only known lysosomal proteins
df = go_tbl_nonan.loc[go_tbl_nonan['KnownLysosomalProteins'] == True].copy()

save_to_folder = r'X:\X\X\X\X'
os.makedirs(save_to_folder, exist_ok=True)

comparisons = {
    '7N/10N': {
        'fc_col': '7N/10N_log2FC',
        'p_col': '7N/10N_pvalue',
        'filename': '7N10N'
    },
    '7T/10T': {
        'fc_col': '7T/10T_log2FC',
        'p_col': '7T/10T_pvalue',
        'filename': '7T10T'
    },
    '7C/10C': {
        'fc_col': '7C/10C_log2FC',
        'p_col': '7C/10C_pvalue',
        'filename': '7C10C'
    }
}

for title, cols in comparisons.items():

    fc_col = cols['fc_col']
    p_col = cols['p_col']
    figurename = cols['filename']

    plot_df = df.dropna(subset=[fc_col, p_col]).copy()

    down = plot_df[
        (plot_df[fc_col] <= -1) &
        (plot_df[p_col] < 0.05)
    ]

    up = plot_df[
        (plot_df[fc_col] >= 1) &
        (plot_df[p_col] < 0.05)
    ]

    plt.figure(figsize=(4.1, 4.8))

    plt.scatter(
        x=plot_df[fc_col],
        y=-np.log10(plot_df[p_col]),
        s=5,
        label="Not significant",
        color="grey"
    )

    plt.scatter(
        x=down[fc_col],
        y=-np.log10(down[p_col]),
        s=5,
        label="Down-regulated",
        color="#1c9099"
    )

    plt.scatter(
        x=up[fc_col],
        y=-np.log10(up[p_col]),
        s=5,
        label="Up-regulated",
        color="#dd1c77"
    )

    plt.xlabel("log2FC")
    plt.ylabel("-log10(p-value)")
    plt.axvline(-1, color="grey", linestyle="--")
    plt.axvline(1, color="grey", linestyle="--")
    plt.axhline(-np.log10(0.05), color="grey", linestyle="--")
    plt.title(title)
    plt.legend()
    plt.tight_layout()

    for extension in ['png', 'eps']:
        plt.savefig(
            os.path.join(save_to_folder, f'{figurename}.{extension}'),
            dpi=300,
            bbox_inches='tight'
        )

    plt.show()