# ##############
#
# # Python script for calculating some properties 
# # reagarding polynya and convection each year
#
# ##############

from myfunctions import *
import gc
from IcePlot import *

from adjustText import adjust_text

def create_new_marker_dict(marker_list, types):
    marker_dict = {}
    if len(marker_list) >= len(types):
        print('create new marker dict')
        for i in range(len(types)):
            marker_dict[types[i]] = marker_list[i]
    else:
        raise Exception("marker list too short")
    return marker_dict


def plot_corr(df_plot, marker_dict, color_dict, save_dir):
    fig, ax = plt.subplots()
    # ax.axvline(c='grey', lw=1)
    # ax.axhline(c='grey', lw=1)
    cnames = df_plot.columns
    cc = 0
    for mtype in df_plot['type'].unique():
        plotdata_df = df_plot.loc[df_plot['type'] == mtype]
    
        ax.plot(plotdata_df[cnames[2]], plotdata_df[cnames[3]], 
                marker=marker_dict[mtype],
                markeredgecolor = color_dict[mtype],
                markerfacecolor = lighten_color(color_dict[mtype], 1.5),
                linestyle='', ms=6, label=mtype)
        cc += 1
    
    texts=[]
    newx=[]
    newy=[]
    for ind in df_plot.index:
        texts+=[ax.text(df_plot[cnames[2]][ind], df_plot[cnames[3]][ind], df_plot['name'][ind], fontsize = 6)]
        newx.append(df_plot[cnames[2]][ind])
        newy.append(df_plot[cnames[3]][ind])

    adjust_text(texts, newx, newy, ax=ax, #avoid_self=True, #force_explode = (0.2, 0.2),
                time_lim=1, #Give it 1 second instead of 0.1 sec to arrange such dense labels
                arrowprops=dict(arrowstyle='-', color='gray', alpha=.5))

    ax.grid(which='both', axis='both',alpha=.5) # draw grid
    ax.set_xlabel(cnames[2])
    ax.set_ylabel(cnames[3])
    ax.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
    savename = save_dir + cnames[2] + '_vs_' + cnames[3] + '.png'
    plt.tight_layout()
    fig.savefig(savename, dpi = 150)


def main():
    import warnings
    warnings.filterwarnings("ignore")
    
    pro1_df = pd.read_csv('properties_diff.csv')
    pro1_df = pro1_df.replace(0.0, np.nan)

    pro2_df = pd.read_csv('properties.csv')
    pro2_df = pro2_df.replace(0.0, np.nan)

    pro_df = pd.merge(pro2_df[pro2_df.columns[0:10]], pro1_df, on="name")

    
    datapd = pd.read_csv('List_model.csv')
    typedf = datapd[['source_id', 'type']]
    typedf = typedf.rename(columns={'source_id': 'name'})
    
    colors = ['#8dd3c7','#ffffb3','#bebada','#fb8072','#d9d9d9','#fdb462','#fccde5','#b3de69','#80b1d3','#bc80bd','#ccebc5']
    marker_list = ['o','v','^','P','X','D','h','*','s','>','<']
    color_dict = create_new_color_dict(colors, datapd, 'type')
    marker_dict = create_new_marker_dict(marker_list, typedf['type'].unique())
    
    for i in range(1, len(pro_df.columns)):
        for j in range(1, len(pro_df.columns)):
            if i < j:
                dfs = pro_df[[pro_df.columns[0], pro_df.columns[i], pro_df.columns[j]]]
                newdf = pd.merge(typedf, dfs, on="name")
                plt_df = newdf.dropna()
                plot_corr(plt_df, marker_dict, color_dict, "Fig_corr_diff/")
                print(i,j)


if __name__ == "__main__":
    main()

