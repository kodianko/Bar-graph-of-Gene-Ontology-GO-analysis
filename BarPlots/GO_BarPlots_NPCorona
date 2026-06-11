import pandas as pd
import matplotlib.pyplot as plt

fbs_path = r'X:\X\X\Programming_death_subroutines_2023_ProteinCoronaData.xlsx'
fbs_tbl = pd.read_excel(fbs_path)


fbs_norm = fbs_tbl[['Normalized FBS', 'Normalized TMA NPs', 'Normalized TMA:MUP MCNPs', 'Normalized TMA:MUA MCNPs', 'Normalized TMA:MUS MCNPs']]

interesting_column_names = ['Albumin', 'AHSG', 'C3', 'Peptidase inhibitor activity', 'Peptidase activity', 'Lipid binding', 'Extracellular matrix',
             'Integrin binding+Actin filament', 'Others']

def col_sum(cols_name):
    return fbs_norm.loc[fbs_tbl[cols_name]].sum()

def barplot_for_one_column(column_name='Albumin', do_plot=True, save_to_folder=False,
                           maximum_height_in_legend=40):
    apo_df = fbs_tbl[['Gene Name', 'Normalized FBS', 'Normalized TMA NPs', 'Normalized TMA:MUP MCNPs', 'Normalized TMA:MUA MCNPs', 'Normalized TMA:MUS MCNPs']].loc[fbs_tbl[
        column_name]]
    apo_df.rename(columns={'Normalized FBS':'FBS', 'Normalized TMA NPs':'TMA', 'Normalized TMA:MUP MCNPs':'TMA:MUP', 'Normalized TMA:MUA MCNPs':'TMA:MUA', 'Normalized TMA:MUS MCNPs':'TMA:MUS'}, inplace = True)
    apo_df.set_index('Gene Name', inplace=True)
    apo_df_T = apo_df.T
    ax = apo_df_T.plot(kind='bar', stacked=True, width=0.7, figsize=(16, 8))
    number_of_proteins_in_legend = apo_df.shape[0]
    number_of_columns_in_legend = 1 + number_of_proteins_in_legend // maximum_height_in_legend
    leg = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., ncol=number_of_columns_in_legend)
    ax.set_position((0.125, 0.15, 0.3, 0.8))
    plt.title(column_name, size = 14)
    plt.ylabel('Normalized abundance')

    # Show sum on each stacked bar
    for i, x in enumerate(col_sum(column_name)):
        ax.annotate(f"{100*x:.01f} %",
                xy=(i, x), xycoords='data',
                xytext=(0, 0.3), textcoords='offset points', ha='center'
                )
    if save_to_folder:
        for extension in ['png', 'eps']:
            ax.get_figure().savefig(f'{save_to_folder}/{column_name}.{extension}', dpi=300)
    if do_plot:
        plt.show()


def barplot_for_many_cols(list_of_col_names, save_to_folder=False):
    for col_name in list_of_col_names:
        barplot_for_one_column(col_name, do_plot=False, save_to_folder=save_to_folder)
    plt.show()

barplot_for_many_cols(interesting_column_names,
                      save_to_folder=r'X:\X\X')

pia_sum = fbs_norm.loc[fbs_tbl['Peptidase inhibitor activity']].sum()
pa_sum = fbs_norm.loc[fbs_tbl['Peptidase activity']].sum()
lb_sum = fbs_norm.loc[fbs_tbl['Lipid binding']].sum()
em_sum = fbs_norm.loc[fbs_tbl['Extracellular matrix']].sum()
intact_sum = fbs_norm.loc[fbs_tbl['Integrin binding+Actin filament']].sum()
others_sum = fbs_norm.loc[fbs_tbl['Others']].sum()
albumin_sum = fbs_norm.loc[fbs_tbl['Albumin']].sum()
ahsg_sum = fbs_norm.loc[fbs_tbl['AHSG']].sum()
c_sum = fbs_norm.loc[fbs_tbl['C3']].sum()

# plot all main groups
go_grp = pd.DataFrame([albumin_sum, ahsg_sum, c_sum, pia_sum, pa_sum, lb_sum, em_sum, intact_sum, others_sum],
                      index=['ALB', 'AHSG', 'C3', 'Peptidase inhibitor activity', 'Peptidase activity', 'Lipid binding', 'Extracellular matrix', 'Integrin binding+Actin filament', 'Others'])
go_grp.rename(columns={'Normalized FBS':'FBS', 'Normalized TMA NPs':'TMA', 'Normalized TMA:MUP MCNPs':'TMA:MUP', 'Normalized TMA:MUA MCNPs':'TMA:MUA', 'Normalized TMA:MUS MCNPs':'TMA:MUS'}, inplace = True)

#go_grp_T = go_grp[fbs_norm].T
go_grp_T = go_grp.T
ax = go_grp_T.plot(kind='bar', stacked=True, width=0.7)
leg = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,)
ax.set_position((0.125,0.15, 0.5, 0.8))
plt.title('all GO groups', size = 14)
plt.show()

save_to_folder=r'X:\X\X'
figurename = "All_GO_groups"
for extension in ['png', 'eps']:
    ax.get_figure().savefig(f'{save_to_folder}/{figurename}.{extension}', dpi=300)

print('gg')
