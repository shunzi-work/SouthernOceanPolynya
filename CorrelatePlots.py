# ##############
#
# # Python script for calculating some properties 
# # reagarding polynya and convection each year
#
# ##############

from myfunctions import *

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from PlotFunctions import create_new_color_dict, create_new_marker_dict, plot_corr, plot_by_type, add_text_annotation, modify_corr_ax, lighten_color
from mycolors import PALETTE as color_palette
from mycolors import MARKERS as marker_list


def plot_corr_from_columns(df_pro, i, j, typedf, marker_dict, color_dict, save_dir, annotext='name'):
    dfs = df_pro[[df_pro.columns[0], df_pro.columns[i], df_pro.columns[j]]]
    newdf = pd.merge(typedf, dfs, on="name")
    plt_df = newdf.dropna()
    plot_corr(plt_df, marker_dict, color_dict, save_dir = save_dir, text = annotext)


def plot_cross_from_df(df_pro, typedf, marker_dict, color_dict, save_dir, annotext='name', onecolumn = None):
    for i in range(1, len(df_pro.columns)):
        for j in range(1, len(df_pro.columns)):
            if i < j:
                if onecolumn:
                    if j != onecolumn:
                        continue
                plot_corr_from_columns(df_pro, i, j, typedf, marker_dict, color_dict, save_dir, annotext=annotext)
                print(i, j)



def plot_figure_p_vs_c(df_pro, df_pro_cross, typedf, marker_dict, color_dict, save_dir, modelname = False):
    dfs = df_pro[['name', 'c_total', 'p_total']]
    dfs2 = df_pro_cross[['name', 'c_cross_rate', 'p_cross_rate']]
    dfs3 = df_pro_cross[['name', 'p_events', 'c_events', 'p_cross_c_events']]

    newdf = pd.merge(typedf, dfs, on="name")
    df_plot = newdf[(newdf['c_total'] != 0) | (newdf['p_total'] != 0)]
    df_plot = df_plot.dropna()

    df_plot2 = pd.merge(typedf, dfs2, on="name")
    df_plot2['p_cross_rate'] = df_plot2['p_cross_rate']* 100
    df_plot2['c_cross_rate'] = df_plot2['c_cross_rate']* 100
    df_plot2 = df_plot2[(df_plot2['c_cross_rate'] != 0) | (df_plot2['p_cross_rate'] != 0)]
    df_plot2 = df_plot2.dropna()

    df_plot3 = pd.merge(typedf, dfs3, on="name")
    df_plot3 = df_plot3[(df_plot3['p_events']!=0)|(df_plot3['c_events']!=0)]
    df_plot3['p_events'] = df_plot3['p_events'] * 100
    df_plot3['c_events'] = df_plot3['c_events'] * 100
    df_plot3['p_cross_c_events'] = df_plot3['p_cross_c_events'] * 100
    df_plot3['new_name'] = df_plot3['number'].astype(str) + ' ' + df_plot3['name']
    df_plot3 = df_plot3.dropna()

    fig = plt.figure(constrained_layout=True, figsize=(7.5, 6.5))
    subfigs = fig.subfigures(nrows=2, ncols=1, height_ratios=[1, 0.8])
    
    ax1, ax2 = subfigs[0].subplots(nrows=1, ncols=2)
    
    plot_by_type(ax1, df_plot, marker_dict, color_dict)
    df11 = df_plot[(df_plot['c_total'] >= 0.3) | (df_plot['p_total']  >= 0.3)]
    add_text_annotation(ax1, df11, t = 'number', fs = 8)
    ax1.set_xlim(-0.5, 10.5)
    modify_corr_ax(ax1, "Total convection region area ($ 10^6\ km^2$)", "Total polynya region area ($ 10^6\ km^2$)", equalbox = True, fs = 10)
    
    left_inset_ax = fig.add_axes([.32, .825, .15, .15], zorder = 1)
    left_inset_ax.set_xlim(-0.025, 0.25)
    left_inset_ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    df12 = df_plot[(df_plot['c_total'] < 0.3) & (df_plot['p_total']  < 0.3)]
    plot_by_type(left_inset_ax, df12, marker_dict, color_dict)
    add_text_annotation(left_inset_ax, df12, t = 'number', fs = 7, expand_rate = 2.5)
    modify_corr_ax(left_inset_ax, '', '', equalbox = True, fs = 8)
    
    
    plot_by_type(ax2, df_plot2, marker_dict, color_dict)
    add_text_annotation(ax2, df_plot2, t = 'number', fs = 8, expand_rate = 2.5)
    ax2.set_xlim(-5, 95)
    modify_corr_ax(ax2, "% of overlap ralative to convection area", "% of overlap ralative to polynya area", equalbox = True, fs = 10)
    
    handles, labels = ax1.get_legend_handles_labels()
    ax2.legend(handles, labels, prop={'size': 9})

    ax3 = subfigs[1].subplots()
    for index, row in df_plot3.iterrows():
        if modelname:
            x = row['new_name']
        else:
            x = row['number'].astype(str)

        
        ax3.bar(
            x, row['p_events'],
            edgecolor = 'k',
            label='polynya' if index == 0 else "", color=color_dict[row['type']], 
            hatch = '//',
            alpha=0.5
        )
        ax3.bar(
            x, row['c_events'], bottom=row['p_events'] - row['p_cross_c_events'], 
            edgecolor = 'k',
            label='convection' if index == 0 else "", color=lighten_color(color_dict[row['type']], 1), alpha=0.5
        )
    ax3.set_ylabel('large event occurance (%)', fontsize=10)
    ax3.grid(which='both', axis='y',alpha=.5) # draw grid
    if modelname:
        ax3.tick_params(axis='x', which='major', labelsize=8, labelrotation=90)
    else:
        ax3.tick_params(axis='x', which='major', labelsize=8)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['bottom'].set_color('#DDDDDD')
    ax3.legend(bbox_to_anchor=(0.65, 0.9), ncol=2, frameon = False, prop={'size': 9})
    # ax3.set_ylim(0, 1.1)
    savename = save_dir + 'p_vs_c.pdf'

    fig.savefig(savename, dpi = 300, bbox_inches='tight')


def main():
    import warnings
    warnings.filterwarnings("ignore")
    
    df_pro = pd.read_csv('properties.csv')

    df_pro_diff = pd.read_csv('properties_diff.csv')
    # df_pro_diff = pd.merge(df_pro[df_pro.columns[0:10]], df_pro_diff, on="name")

    df_pro_cross = pd.read_csv('CrossRate.csv')

    df_pro_fft = pd.read_csv('periodicity.csv')
    df_pro_fft = df_pro_fft.replace(0.0, np.nan)
    
    df_pro = df_pro.replace(0.0, np.nan)
    df_pro_diff = df_pro_diff.replace(0.0, np.nan)

    df_pro_ncdiff = pd.read_csv('properties_pcdiff.csv')
    # df_pro_ncdiff = pd.merge(df_pro[df_pro.columns[0:10]], df_pro_ncdiff, on="name")
    df_pro_ncdiff = df_pro_ncdiff.replace(0.0, np.nan)

    # df_pro_plt_new = pd.merge(df_pro, df_pro_cross[df_pro_cross.columns[0:2]], on="name")
    df_all1 = pd.merge(df_pro, df_pro_diff, on="name")
    df_all2 = pd.merge(df_all1, df_pro_cross, on="name") 
    df_all3 = pd.merge(df_all2, df_pro_fft, on="name")
    df_all = pd.merge(df_all3, df_pro_ncdiff, on="name")

    datapd = pd.read_csv('List_model.csv')
    typedf = datapd[['source_id', 'type']]
    typedf['number'] = datapd.index + 1
    typedf = typedf.rename(columns={'source_id': 'name'})
    typedf['long_name'] = typedf['number'].astype(str) + ' ' + typedf['name']
    
    color_dict = create_new_color_dict(color_palette, typedf['type'].unique())
    marker_dict = create_new_marker_dict(marker_list, typedf['type'].unique())

    save_dir = "../../figure_out/CorrelatePlots2/"

    # print("Plotting correlations for properties...")
    # plot_cross_from_df(df_pro, typedf, marker_dict, color_dict, save_dir, annotext='number')

    # print("Plotting correlations for property differences...")
    # plot_cross_from_df(df_pro_diff, typedf, marker_dict, color_dict, save_dir, annotext='number')

    # print("Plotting correlations for property differences in n/c/p...")
    # plot_cross_from_df(df_pro_ncdiff, typedf, marker_dict, color_dict, save_dir, annotext='number')

    # print("Plotting correlations for periodicity...")
    # plot_corr_from_columns(df_pro_fft, 1, 3, typedf, marker_dict, color_dict, save_dir, annotext='number')

    # print("Plotting correlations for cross rates...")
    # plot_cross_from_df(df_pro_cross, typedf, marker_dict, color_dict, save_dir, annotext='number')

    print("Plotting correlations...")
    plot_cross_from_df(df_all, typedf, marker_dict, color_dict, save_dir, annotext='number')

    # plot_figure_p_vs_c(df_pro, df_pro_cross, typedf, marker_dict, color_dict, save_dir, modelname = False)


if __name__ == "__main__":
    main()

