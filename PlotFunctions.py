###############
# 
# # Python script for plotting
# 
###############

from myfunctions import *
import numpy as np
import gc

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.colors as mcolors
import colorsys

import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

import cmocean

from CalculateProperties import open_mld
from adjustText import adjust_text

def calculate_mean_ice(ice, area):
    return (ice*area).sum()/area.where(ice>=0).sum()

def calculate_pmean_ice(ice, area):
    icemax = ice.max("time")
    mean = (ice.where(icemax>0)*area).sum()/(area.where(icemax>0)).sum()/len(ice.time)
    return mean

def lighten_color(color, amount=0.5):
    """
    from https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    
    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    c = np.array(colorsys.rgb_to_hls(*mcolors.to_rgb(color)))
    return colorsys.hls_to_rgb(c[0],1-amount * (1-c[1]),c[2])

def plot_ice_hist(datapd, pice, p_count, savename, figsize=(6.5, 7)):
    fig = plt.figure(figsize=figsize)
    n = 0
    ice_thresholds = np.arange(0, 100, step=1)
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        if ispickleexists(name, p_count) and ispickleexists(name, pice):
            n+=1
            ds = np.array(openpickle(name, p_count))
            dssiconc = openpickle(name, pice)
            ice_mean = calculate_mean_ice(dssiconc.siconc, dssiconc.areacello)
            ice_mean_not0 = calculate_pmean_ice(dssiconc.siconc, dssiconc.areacello)
     
            ax = fig.add_subplot(10, 5, n)
            plt.subplots_adjust(left=0.015,
                                bottom=0.055, 
                                right=0.985, 
                                top=0.99, 
                                wspace=0.08, 
                                hspace=0.08)
            
            ax.set_ylim((-0.05, 1.05))
            ax.tick_params(labelsize = 6)
            
            titletext = "{}, {:.1f}".format(name, np.max(ds[1:])/1e12)
            ax.set_title(titletext, fontsize=6, y=1.0, pad=-4)
            # ax.text(5, 0.7, text_mean, fontsize = 6)
            
            ax.set_xlim((-1,101))
            ax.set_xticks([0,20,40,60,80,100])
    
            ax2 = ax.twinx()  # instantiate a second Axes that shares the same x-axis
            icevalues = dssiconc.siconc.values.flatten()
            icev = icevalues[icevalues>0]
            
            pa = ax.plot(ice_thresholds[:], ds[:,0]/np.max(ds[:]), linewidth = 1, 
                         color = 'b', label = 'polynya area')
            # pa2 = ax.plot(ice_thresholds[:], ds[:,1]/np.max(ds[:]), linewidth = 1, 
            #              color = 'r', label = 'polynya area')
            # pa3 = ax.plot(ice_thresholds[:], ds[:,2]/np.max(ds[:]), linewidth = 1, 
            #              color = 'y', label = 'polynya area')
            pc = ax2.axvline(x=ice_mean.values.item(), color = 'g', linestyle = '--',
                            linewidth=1, label = 'mean SIC in SO')
            pe = ax.axvline(x=ice_mean_not0.values.item(), color = 'orange', linestyle = '-.',
                            linewidth=1, label = 'mean SIC within SIE')
            pb = ax2.hist(icev, bins=100, color = 'red', edgecolor=None,
                          alpha = 0.6, label = 'sea ice concentration')
    
            if n<=44:
                ax.set_xticklabels([])            
            if n == 48:
                ax.set_xlabel('ice concentration (%)', fontsize=10)
                # ax.xaxis.set_label_coords(1.6, -0.35)
                
            ax.set_yticklabels([])
            ax.set_frame_on(False)
            ax.tick_params(tick1On=False)
            ax2.set_yticklabels([])
            ax2.set_yticks([])
            ax2.set_frame_on(False)
            ax2.tick_params(tick1On=False)
            
            ax.xaxis.grid(True,'major', ls='-', lw=0.2, alpha = 0.6)
            if n == 49:
                l = ax.legend(fontsize = 6, frameon=False, loc='center left', bbox_to_anchor=(1.05, 0.4))
                l2 = ax2.legend(fontsize = 6, frameon=False, loc='center left', bbox_to_anchor=(1.05, 0))
            # break
    fig.savefig(savename, format='pdf')


def create_new_color_dict(color_palette, types, lighten = 1):
    color_dict = {}
    if len(color_palette) >= len(types):
        print('create new color dict')
        for i in range(len(types)):
            color_dict[types[i]] = lighten_color(color_palette[i], lighten)
        return color_dict
    else:
        raise ValueError("More types than colors.")
    
def create_new_marker_dict(marker_list, types):
    marker_dict = {}
    if len(marker_list) >= len(types):
        print('create new marker dict')
        for i in range(len(types)):
            marker_dict[types[i]] = marker_list[i]
        return marker_dict
    else:
        raise ValueError("marker list too short")

def modify_map(ax):
    circle = create_circle()
    ax.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, zorder=1, color = 'grey')
    # ax.add_feature(cfeature.COASTLINE, linewidth = 0.5)
    ax.add_feature(cfeature.OCEAN, alpha = 0.15)
    ax.set_boundary(circle, transform=ax.transAxes)
    # ax.spines['geo'].set_linewidth(0.5)
    ax.spines['geo'].set_edgecolor(None)
    # ax.gridlines(draw_labels=False, 
    #              ylocs=np.linspace(-90, 90, 7), 
    #              color = 'grey', linestyle = '-.', linewidth = 0.5, alpha = 0.8)


def plot_type_color(fig, ax, color, titleonly):
    bbox = ax.get_position()
    if titleonly:
        rect = Rectangle((bbox.x0-0.002,bbox.y1+0.001),
                         bbox.width+0.004, 0.012, 
                         fill=True, color=color, alpha=1, zorder=-1,
                         transform=fig.transFigure, clip_on=False)
    else:
        rect = Rectangle((bbox.x0-0.002367,bbox.y0),
                         bbox.width+0.004734,bbox.height+0.02, 
                         fill=True, color=color, alpha=1, zorder=-1,
                         transform=fig.transFigure, clip_on=False)
    fig.add_artist(rect)

def create_circle():
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    return circle

def get_icemax_polynya(name, p_ice, ice_threshold):
    dsice = openpickle(name, p_ice)
    plt_icemax = dsice.siconc.max('time')
    pltx = dsice.newlon
    plty = dsice.newlat
    flood_points = [(0,0)]
    if name == "MRI-ESM2-0":
        flood_points = [(0,0), (0,40), (0,100), (0,200)]
    mask = detect_polynya(dsice.siconc, dsice.areacello, ice_threshold, flood_points = flood_points)
    plt_polynya = mask.count('time').where(mask.count('time')>0)/len(mask.time)
    return pltx, plty, plt_icemax, plt_polynya, len(dsice.time)

def add_subplot_icepolynya(fig, n, name, color, titleonly, reso_text, pltx,plty,plt_icemax,plt_polynya, timelength):
    ax = fig.add_subplot(7, 8, n, projection=ccrs.SouthPolarStereo())
    plt.subplots_adjust(left=0.01,
                        bottom=0.01, 
                        right=0.99, 
                        top=0.99, 
                        wspace=0.04, 
                        hspace=0.04)
    plot_type_color(fig, ax, color, titleonly)
    modify_map(ax)
    add_subplot_text(ax, name, n, timelength, reso_text)

    im = ax.pcolormesh(pltx, plty, plt_icemax, vmin = 0, vmax=100, 
                        transform=ccrs.PlateCarree(), cmap=cmocean.cm.ice)
    im2 = ax.pcolormesh(pltx, plty, plt_polynya, vmin = 0, vmax=0.2, 
                        transform=ccrs.PlateCarree(), cmap=plt.cm.plasma)
    return im, im2

def add_subplot_mld(fig, n, name, color, titleonly, reso_text, pltx, plty, plt_mld, timelength):
    ax = fig.add_subplot(7, 8, n, projection=ccrs.SouthPolarStereo())
    plt.subplots_adjust(left=0.01,
                        bottom=0.01, 
                        right=0.99, 
                        top=0.99, 
                        wspace=0.04, 
                        hspace=0.04)
    modify_map(ax)
    plot_type_color(fig, ax, color, titleonly)
    add_subplot_text(ax, name, n, timelength, reso_text)
    im = ax.pcolormesh(pltx, plty, plt_mld, vmin = 0, vmax=0.2, 
                        transform=ccrs.PlateCarree(), cmap=plt.cm.plasma)
    return im


def to_percentage(x, pos):
    return f"{x * 100:.0f}%"

def add_cbar(fig, axes_loc, im, label_text, format_percent=False):
    cbar_ax1 = fig.add_axes(axes_loc)
    cbar1 = fig.colorbar(im, cax=cbar_ax1, orientation='horizontal')
    cbar1.set_label(label_text, size=8, labelpad=-0.1)
    cbar1.ax.tick_params(labelsize=6, direction='in')
    cbar1.outline.set_visible(False)
    if format_percent:
        cbar1.formatter = FuncFormatter(to_percentage)
        cbar1.update_ticks()

def add_cbars(fig, im, im2):
    add_cbar(fig, [0.62, 0.05, 0.35, 0.01], im, 'Sea ice concentration (%)')
    add_cbar(fig, [0.62, 0.11, 0.35, 0.01], im2, 'Frequency of occurrence', format_percent = True)

def add_type_color_legend(fig, color_dict, legendname, bta = (0.55, 0.135), n_cols = 3, fs = 6):
    types = list(color_dict.keys())
    proxy_artists = [Rectangle((0,0),1,1, color = color_dict[c]) for c in types]
    fig.legend(proxy_artists, types, 
               title = legendname,#'ocean module types',
               frameon = False,
               bbox_to_anchor=bta,
               ncols=n_cols, fontsize = fs)

def plot_polynya_maps(
    datapd, p_ice, color_palette, ice_threshold, 
    typename='type_ice', titleonly=True, 
    resotext=False, lighten =1, figsize=(6.5,7.5),
    savename_pre = 'Figures/polynya_maps_'
):
    fig = plt.figure(figsize=figsize)
    n = 0
    color_dict = create_new_color_dict(color_palette, datapd[typename].unique(), lighten = lighten)
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        t = datapd.at[i, typename]
        pltx, plty, plt_icemax, plt_polynya, timelength = get_icemax_polynya(name, p_ice, ice_threshold)
        n = n + 1
        if resotext:
            reso_text = str(datapd.at[i, 'resolution'])
        else:
            reso_text = None
        im, im2 = add_subplot_icepolynya(fig, n, name, color_dict[t], titleonly, reso_text, pltx,plty,plt_icemax,plt_polynya, timelength)
        gc.collect()
    add_cbars(fig, im, im2)
    add_type_color_legend(fig, color_dict, 'sea ice module types')
    text_threshold = str(ice_threshold) + '%'
    fig.text(0.15,0.06, text_threshold)
    figsavename = savename_pre + str(ice_threshold) + '.png'
    fig.savefig(figsavename, dpi = 300)


def add_subplot_text(ax, name, n, timelength, resotext):
    title = ax.set_title('{}'.format(name), 
                             fontsize=6, pad=-0.5)
    # title._bbox_patch._mutation_aspect = 1
    # title.get_bbox_patch().set_boxstyle("square", pad=11.9)
    ax.text(0,-90, timelength, fontsize=6, color='w', ha='center')
    ax.text(180,-55, resotext, transform=ccrs.PlateCarree(), fontsize=6, color='k', ha='center')
    ax.text(45,-47, n, transform=ccrs.PlateCarree(), fontsize=6, color='k', ha='center')


def plot_polynya_maps_from_polynya_data(datapd, p_ice, p_polynya, colors, typename='type_ice', titleonly=True, resotext=False, lighten =1, figsize=(6.5,7)):
    fig = plt.figure(figsize=figsize)
    n = 0
    color_dict = create_new_color_dict(colors, datapd[typename].unique(), lighten = lighten)
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        t = datapd.at[i, typename]
        dsice = openpickle(name, p_ice)
        dspolynya = openpickle(name, p_polynya)
        pltx, plty = dsice.newlon, dsice.newlat

        plt_icemax = dsice.siconc.max('time')
        polynya_count = dspolynya.where(dspolynya>0).count('time')
        timelength = len(dspolynya.time)
        plt_polynya = polynya_count.where(polynya_count>0)/timelength

        n = n + 1
        if resotext:
            reso_text = str(datapd.at[i, 'resolution'])
        else:
            reso_text = None
        im, im2 = add_subplot_icepolynya(fig, n, name, color_dict[t], titleonly, reso_text, pltx,plty,plt_icemax,plt_polynya, timelength)
        gc.collect()
    add_cbars(fig, im, im2)
    add_type_color_legend(fig, color_dict, 'sea ice module types')
    
    figsavename = 'polynya_maps.png'
    fig.savefig(figsavename, dpi = 300)


def plot_convection_maps(datapd, p_mld, p_mlotst, colors, typename='type', titleonly=True, resotext=True, lighten =1, figsize=(6.5,7)):
    fig = plt.figure(figsize=figsize)
    n = 0
    color_dict = create_new_color_dict(colors, datapd[typename].unique(), lighten = lighten)
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        t = datapd.at[i, typename]
        damld, dsmld = open_mld(p_mlotst, p_mld, name)
        if len(damld.time)>500:
            damld = damld.isel(time = slice(-500, None))
        
        pltx = dsmld.newlon
        plty = dsmld.newlat
        plt_mld_count = damld.where(damld>=2000).count("time")
        timelength = len(damld.time)
        plt_mld = plt_mld_count.where(plt_mld_count>0)/timelength

        n = n + 1
        if resotext:
            reso_text = str(datapd.at[i, 'resolution'])
        else:
            reso_text = None
        im = add_subplot_mld(fig, n, name, color_dict[t], titleonly, reso_text, pltx, plty, plt_mld, timelength)
        gc.collect()
    add_cbar(fig, [0.62, 0.08, 0.35, 0.01], im, 'Frequency of occurrence', format_percent = True)
    add_type_color_legend(fig, color_dict, 'ocean module types')
    fig.savefig("convection_maps.png", dpi = 300)


def plot_by_type(ax, df_plot, marker_dict, color_dict, ms = 6, col1 = -2, col2 = -1):
    var1 = df_plot.columns[col1]
    var2 = df_plot.columns[col2]
    for mtype in df_plot['type'].unique():
        plotdata_df = df_plot.loc[df_plot['type'] == mtype]
        ax.plot(plotdata_df[var1], plotdata_df[var2], 
                marker = marker_dict[mtype],
                markerfacecolor = lighten_color(color_dict[mtype], 1.5),
                markeredgecolor = lighten_color(color_dict[mtype], 2),
                linestyle = '',
                ms=ms, label=mtype)


def add_text_annotation(ax, df_plot, t = 'name', fs = 6, expand_rate = 2, col1 = -2, col2 = -1):
    var1 = df_plot.columns[col1]
    var2 = df_plot.columns[col2]
    
    texts=[]
    newx=[]
    newy=[]
    for ind in df_plot.index:
        texts+=[ax.text(df_plot[var1][ind], df_plot[var2][ind], df_plot[t][ind], fontsize = fs, alpha = 0.6)]
        newx.append(df_plot[var1][ind])
        newy.append(df_plot[var2][ind])

    adjust_text(texts, newx, newy, ax=ax, avoid_self=True, 
                #force_explode = (0.2, 0.2),
                expand=(expand_rate, expand_rate),
                time_lim = 1, #Give it 1 second instead of 0.1 sec to arrange such dense labels
                arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5))
    

def modify_corr_ax(ax, xlabel, ylabel, equalbox = False, fs = 8):
    ax.grid(which='both', axis='both',alpha=.5) # draw grid
    ax.set_xlabel(xlabel, fontsize = fs)
    ax.set_ylabel(ylabel, fontsize = fs)
    if equalbox:
        ax.set_aspect('equal')
        ax.set_box_aspect(1)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    ax.tick_params(axis='both', which='minor', labelsize=fs)
    

def plot_corr(
    df_plot, 
    marker_dict, 
    color_dict, 
    vline = False, 
    hline = False, 
    text = 'name',
    equalbox = False, 
    save_dir = 'Fig_corr/'
):
    fig, ax = plt.subplots()
    if vline:
        ax.axvline(c='grey', lw=1)
    if hline:
        ax.axhline(c='grey', lw=1)

    plot_by_type(ax, df_plot, marker_dict, color_dict)
    add_text_annotation(ax, df_plot, t = text)

    modify_corr_ax(ax, xlabel=df_plot.columns[-2], ylabel=df_plot.columns[-1], equalbox = equalbox)

    ax.legend(bbox_to_anchor=(1.04, 0.5), loc="center left")
    savename = save_dir + df_plot.columns[-2] + '_vs_' + df_plot.columns[-1] + '.png'
    plt.tight_layout()
    fig.savefig(savename, dpi = 150)