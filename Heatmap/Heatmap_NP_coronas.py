import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy

# =========================
# Load data
# =========================
df1_path = r'X:\X\X\X\NP_coronas_data.xlsx'
df1_tbl = pd.read_excel(df1_path)

fer_path = r'X:\X\X\X\NP_coronas_proteins_groups_for_heatmap.xlsx'
df2_tbl = pd.read_excel(fer_path)

# =========================
# Clean gene names
# =========================
df1_tbl['Gene_clean'] = df1_tbl['Gene Name'].astype(str).str.upper().str.strip()
df2_tbl['Gene_clean'] = df2_tbl['Gene'].astype(str).str.upper().str.strip()

# Remove duplicates in df2 (important)
df2_tbl = df2_tbl.drop_duplicates('Gene_clean')

# =========================
# Create ordering from df2
# =========================
order_dict = {gene: i for i, gene in enumerate(df2_tbl['Gene_clean'])}

# =========================
# Merge (bring categories)
# =========================
df_merged = df1_tbl.merge(
    df2_tbl[['Gene_clean', 'Protein Class']],
    how='left',
    on='Gene_clean'
)

# =========================
# Apply df2 ordering
# =========================
df_merged['order_from_df2'] = df_merged['Gene_clean'].map(order_dict)

df_merged = df_merged.sort_values(
    by='order_from_df2',
    na_position='last'
)

# =========================
# Add Ferroptosis column
# =========================
df_merged['Protein group'] = df_merged['Protein Class'].notna()

# =========================
# Preserve category order (CRITICAL)
# =========================
cat_order = pd.unique(df2_tbl['Protein Class'])

df_merged['Protein Class'] = pd.Categorical(
    df_merged['Protein Class'],
    categories=cat_order,
    ordered=True
)

# =========================
# Cleanup
# =========================
df_merged = df_merged.drop(columns=['Gene_clean', 'order_from_df2'])
# df_merged.columns = df_merged.columns.get_level_values(0)
# =========================
# Save result
# =========================
# output_path = r'E:\!!!Lysosomal Proteomics\Lysosomal proteins\with_Ferroptosis_merge.xlsx'
# df_merged.to_excel(output_path, index=False)

non_nan = df_merged.dropna(subset=['Relative1 FBS', 'Relative2 FBS','Relative3 FBS','Relative1 TMA NPs','Relative2 TMA NPs', 'Relative3 TMA NPs',
                                   'Relative1 TMA:MUA MCNPs', 'Relative2 TMA:MUA MCNPs', 'Relative3 TMA:MUA MCNPs', 'Relative1 TMA:MUS MCNPs', 'Relative2 TMA:MUS MCNPs', 'Relative3 TMA:MUS MCNPs', 'Relative1 TMA:MUP MCNPs', 'Relative2 TMA:MUP MCNPs', 'Relative3 TMA:MUP MCNPs'])

np_tbl_topnan = df_merged.loc[
    df_merged['Protein group'] &
    (df_merged['Accession'].astype(str).str.strip() != 'P02070')
]

sample_cols = ['FBS1', 'FBS2', 'FBS3', 'TMA1','TMA2', 'TMA3', 'MUA1','MUA2', 'MUA3', 'MUS1', 'MUS2', 'MUS3', 'MUP1','MUP2', 'MUP3']
col = np_tbl_topnan.columns[list(range(37,40)) + list(range(9,12)) + list(range(23,26)) + list(range(30,33)) + list(range(16,19))]

matrix = np_tbl_topnan[col].copy()
ytick = np_tbl_topnan['Gene Name']
labels = np_tbl_topnan['Protein Class']

lut = dict(zip(set(labels), sns.hls_palette(len(set(labels)), l=0.5, s=0.8)))
labels_clean = labels.astype(str)
row_colors = labels_clean.map(lut)

prot_class = labels_clean.to_list()
prot_class_getind = sorted(np.unique(prot_class, return_index=True)[1])
prot_class_names = pd.unique(np.array(prot_class)).tolist()

g = sns.clustermap(matrix, col_cluster=False, row_cluster=False,
                 cmap='coolwarm',
                 row_colors=row_colors,
                 yticklabels=ytick,
                 xticklabels=sample_cols,
                 norm=LogNorm(vmin=matrix.min().min(), vmax=matrix.max().max()),
                 figsize=(10, 12),
                 cbar_pos=[.28, .98, .625, .02],  # x, y, width, height in "figure coordinates"
                 cbar_kws={'orientation': "horizontal", 'label': 'Normalized Abundance ($log_{10}$)'}
                 )

hm = g.ax_heatmap.get_position()
g.ax_heatmap.set_position([hm.x0*1.33, hm.y0*0.90, hm.width*0.9, hm.height*1.18])

g.ax_row_colors.set_yticks(0.5 * (np.array(prot_class_getind) + np.array(prot_class_getind[1:] + [len(prot_class)])))
g.ax_row_colors.set_yticklabels(prot_class_names, fontsize=12)
g.ax_row_colors.yaxis.set_tick_params(size=0) # make tick marks invisible

class_bar = g.ax_row_colors.get_position()
g.ax_row_colors.set_position([class_bar.x0*1.42, class_bar.y0*0.90, class_bar.width*0.8, class_bar.height*1.18])
g.ax_row_colors.set_title(label=None)

# plt.show()

save_to_folder = r'X:\X\X\V'
figurename = "Top50"
for extension in ['png', 'eps']:
    g.fig.savefig(f'{save_to_folder}/{figurename}.{extension}', dpi=300)

plt.show()
