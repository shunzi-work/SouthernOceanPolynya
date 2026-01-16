from myfunctions import *
from PlotFunctions import *

from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

def detect_polynya_sithick(da_ice, da_thick, da_area, ice_threshold, area_threshold=(100, 1000), flood_points = [(0,0)], buffering = 15):
    from scipy import ndimage
    from skimage.segmentation import flood_fill
    # Copy the input data arrays to avoid modifying the original data
    daice = da_ice.copy()
    dathick = da_thick.copy()
    daarea = da_area.copy()

    # Fill NaN values in the area data array with the mean value
    daarea = daarea.fillna(daarea.mean().values.item())
    
    # Generate a binary structure for labeling connected components
    s = ndimage.generate_binary_structure(2,2)
    
    # Create an empty masked data array with the same dimensions as the input data
    da_masked = xr.DataArray(np.nan*np.empty_like(dathick), dims = dathick.dims, coords = dathick.coords)
    
    # Loop through each time step in the data array
    for year in daice.time:
        # Copy the sea ice concentration data for the current time step and fill NaN values with 0
        ice0_flood = daice.sel(time = year).copy()
        ice0_flood = ice0_flood.fillna(0).values
        
        # Apply flood fill to remove coastal areas based on the provided flood points
        for flood_point in flood_points:
            ice0_flood = flood_fill(ice0_flood, flood_point, 0, tolerance=buffering)
        
        # Apply flood fill to remove areas outside the sea ice extent
        ice0_flood = flood_fill(ice0_flood, (ice0_flood.shape[0]-1, ice0_flood.shape[1]-1), 0, tolerance=buffering)

        icethick0 = dathick.sel(time = year).copy()
        icethick0 = icethick0.fillna(0).values
        icethick0 = np.where(ice0_flood>0, icethick0, 0)

        # Identify areas with sea ice concentration below the threshold
        icenew = icethick0 <= ice_threshold
        
        # Label connected components in the binary image
        labeled_image, num_features = ndimage.label(icenew, structure = s)
        
        # If there are less than 2 features, skip to the next time step
        if num_features < 2:
            continue
        
        # Initialize a mask to store the detected polynyas
        mask = np.zeros_like(labeled_image)
        
        # Loop through each labeled feature
        for i in range(1, num_features+1):
            # Calculate the area of the current feature
            area = daarea.where(labeled_image == i).sum()/1e9  # m2 -> 10^3 km2
            
            # If the area is within the specified thresholds, mark it as a polynya
            if (area > area_threshold[0]) and (area < area_threshold[1]):
                mask[labeled_image == i] = 1
        
        # Mask the sea ice concentration data for the current time step
        ice_value = dathick.sel(time=year).values
        ice_value[mask == 0] = np.nan
        da_masked.loc[year] = ice_value
    
    # Return the masked data array
    return da_masked


def adjust_subplot():
    plt.subplots_adjust(left=0.015,
                        bottom=0.02, 
                        right=0.925,
                        top=0.98,
                        wspace=0.15,
                        hspace=0.08)


def plot_methods(data_mld, data_siconc, data_sithick, figsize=(4.5, 3.5), savename = 'Methods.png'):
    fig = plt.figure(figsize=figsize)
    
    ax1 = fig.add_subplot(2, 2, 1, projection=ccrs.SouthPolarStereo())
    adjust_subplot()
    
    mldavg = data_mld.mld.mean('time')
    im1 = ax1.pcolormesh(data_mld.newlon, data_mld.newlat, mldavg, 
                      transform=ccrs.PlateCarree(), cmap=cmocean.cm.deep)
    im1.set_clim(0, 1000)
    cbar1 = plt.colorbar(im1, ax=ax1, label='meters', pad=0.05)
    cbar1.set_label('Mean Sept. MLD (m)', size=8)
    cbar1.ax.tick_params(labelsize=6)
    cbar1.outline.set_visible(False)
    
    mld1000 = xr.where(data_mld.mld.where(data_mld.mld>=1000).mean('time')>0, 1, np.nan)
    ax1.pcolormesh(data_mld.newlon, data_mld.newlat, mld1000,
                   transform=ccrs.PlateCarree(), cmap=mcolors.ListedColormap(['orange']), alpha=0.5)
    
    mld2000 = xr.where(data_mld.mld.where(data_mld.mld>=2000).mean('time')>0, 1, np.nan)
    ax1.pcolormesh(data_mld.newlon, data_mld.newlat, mld2000,
                   transform=ccrs.PlateCarree(), cmap=mcolors.ListedColormap(['red']), alpha=0.3)
    ax1.text(-45, -40, 'a.', transform=ccrs.PlateCarree(), fontweight='bold')
    modify_map(ax1)
    
    
    ax2 = fig.add_subplot(2, 2, 2, projection=ccrs.SouthPolarStereo())
    adjust_subplot()
    mldstd = data_mld.mld.std("time")
    im2 = ax2.pcolormesh(data_mld.newlon, data_mld.newlat, mldstd, 
                         transform=ccrs.PlateCarree(), cmap=cmocean.cm.amp)
    cbar2 = plt.colorbar(im2, ax=ax2, label='meters', pad=0.05)
    cbar2.set_label('Sept. MLD STD (m)', size=8)
    cbar2.ax.tick_params(labelsize=6)
    cbar2.outline.set_visible(False)
    mldzscore = (mldstd - mldstd.mean())/mldstd.std()
    ax2.contour(data_mld.newlon, data_mld.newlat, mldzscore, 
                transform=ccrs.PlateCarree(),
                levels=[2], colors='k', linewidths=0.5)
    ax2.text(-45, -40, 'b.', transform=ccrs.PlateCarree(), fontweight='bold')
    modify_map(ax2)
    
    
    ax3 = fig.add_subplot(2, 2, 3, projection=ccrs.SouthPolarStereo())
    adjust_subplot()
    siconcmax = data_siconc.siconc.max('time')
    im3 = ax3.pcolormesh(data_siconc.newlon, data_siconc.newlat, siconcmax,
                         transform=ccrs.PlateCarree(), cmap=cmocean.cm.ice)
    
    cbar3 = plt.colorbar(im3, ax=ax3, label='fraction', pad=0.05)
    cbar3.set_label('Max Sept. SIC (%)', size=8)
    cbar3.ax.tick_params(labelsize=6)
    cbar3.outline.set_visible(False)
    
    
    ax3.contour(data_siconc.newlon, data_siconc.newlat, 
                xr.where(detect_polynya(data_siconc.siconc, data_siconc.areacello, 15).count('time')>0, 1.0, 0),
                transform=ccrs.PlateCarree(), colors='C1', linewidths=1, alpha = 0.2)
    
    ax3.contour(data_siconc.newlon, data_siconc.newlat, 
                xr.where(detect_polynya(data_siconc.siconc, data_siconc.areacello, 45).count('time')>0, 1.0, 0),
                transform=ccrs.PlateCarree(), colors='C3', linewidths=1, alpha = 0.2)
    
    ax3.contour(data_siconc.newlon, data_siconc.newlat, 
                xr.where(detect_polynya(data_siconc.siconc, data_siconc.areacello, 75).count('time')>0, 1.0, 0),
                transform=ccrs.PlateCarree(), colors='C9', linewidths=1, alpha = 0.2)
    
    ax3.text(-45, -40, 'c.', transform=ccrs.PlateCarree(), fontweight='bold')
    modify_map(ax3)
    
    
    ax4 = fig.add_subplot(2, 2, 4, projection=ccrs.SouthPolarStereo())
    adjust_subplot()
    sithickavg = data_sithick.sithick.mean('time')
    im4 = ax4.pcolormesh(data_sithick.newlon, data_sithick.newlat, sithickavg.where(sithickavg>0),
                         transform=ccrs.PlateCarree(), cmap=cmocean.cm.rain)
    im4.set_clim(0, 1.5)
    ax4.contour(data_sithick.newlon, data_sithick.newlat, 
                xr.where(detect_polynya_sithick(data_siconc.siconc, data_sithick.sithick, data_sithick.areacello, 0.2).count('time')>0, 1.0, 0),
                transform=ccrs.PlateCarree(), colors='C1', linewidths=1, alpha = 0.2)
    ax4.contour(data_sithick.newlon, data_sithick.newlat, 
                xr.where(detect_polynya_sithick(data_siconc.siconc, data_sithick.sithick, data_sithick.areacello, 0.5).count('time')>0, 1.0, 0),
                transform=ccrs.PlateCarree(), colors='C3', linewidths=1, alpha = 0.2)
    ax4.contour(data_sithick.newlon, data_sithick.newlat, 
                xr.where(detect_polynya_sithick(data_siconc.siconc, data_sithick.sithick, data_sithick.areacello, 0.8).count('time')>0, 1.0, 0),
                transform=ccrs.PlateCarree(), colors='C9', linewidths=1, alpha = 0.2)
    
    cbar4 = plt.colorbar(im4, ax=ax4, label='meters', pad=0.05)
    cbar4.set_label('Mean Sept.\n Sea Ice Thickness (m)', size=8)
    cbar4.ax.tick_params(labelsize=6)
    cbar4.outline.set_visible(False)
    ax4.text(-45, -40, 'd.', transform=ccrs.PlateCarree(), fontweight='bold')
    modify_map(ax4)
    
    fig.savefig(savename, dpi = 300)
    # fig.savefig('Methods.eps', format='eps')  # too big (should be less than 5MB
    # fig.savefig('Methods.pdf', format='pdf')  # too big (should be less than 5MB


def plot_detection(ds_ice, iyear = 188, fsize = (6, 4.7), savename = 'polynya_detect.png'):
    ice_year = ds_ice.siconc.isel(time = slice(iyear, iyear+1))
    
    fig = plt.figure(constrained_layout=True, figsize=fsize)
    subfigs = fig.subfigures(nrows=2, ncols=1, height_ratios=[0.6, 1])
    
    area_p = []
    for t in range(1, 99):
        p_i = detect_polynya(ice_year, ds_ice.areacello, t, area_threshold=[0, 2000])
        area = ds_ice.areacello.where(p_i.isel(time = 0) > 0).sum()/1e12
        area_p.append(area.item())
    
    gs = GridSpec(1, 6, figure=subfigs[0])
    ax0 = subfigs[0].add_subplot(gs[0, 2:])
    ax0.plot(range(1, 99), area_p, color = 'k')
    
    ax0.grid("on", alpha = 0.3)
    plt.box(False)
    ax0.set_xlabel("SIC threshold (%)")
    ax0.set_ylabel("polynya area\n($10^6\ km^2$)")
    ax0.set_xticks(range(0, 110, 10))
    ax0.set_yticks(np.arange(-0.5, 2.5, 0.5))
    ax0.text(-32, 1.9, "a.", fontweight='bold')
    ax0.text(98, 1.9, "b.", fontweight='bold')
    ax0.text(98, -1.7, "c.", fontweight='bold')
    ax0.set_xlim([0, 100])
    ax0.set_ylim([-0.1, 2.1])
    
    map_ax = fig.add_axes([0.065, 0.70, 0.275, 0.275], zorder=1, projection=ccrs.SouthPolarStereo())
    modify_map(map_ax)
    im = map_ax.pcolormesh(ds_ice.newlon, ds_ice.newlat, ice_year.isel(time = 0), vmin = 0, vmax = 100,
                           transform=ccrs.PlateCarree(), cmap=cmocean.cm.ice)
    cbar_ax1 = fig.add_axes([0.01, 0.70, 0.01, 0.275])
    cbar1 = fig.colorbar(im, cax=cbar_ax1, orientation='vertical')
    cbar1.set_label("SIC (%)", labelpad=-0.1)
    cbar1.ax.tick_params(direction='in')
    cbar1.outline.set_visible(False)
    
    
    for i in range(0, 8):
        test = detect_polynya(ice_year, ds_ice.areacello, 10*i+10, area_threshold=[0, 2000])
        test = xr.where(test > 0, 1, np.nan)
        ax = subfigs[1].add_subplot(2, 4, i+1, projection=ccrs.SouthPolarStereo())
        plt.subplots_adjust(left=0.01,
                                bottom=0.01, 
                                right=0.99, 
                                top=0.99, 
                                wspace=0.04, 
                                hspace=0.04)
        modify_map(ax)
        
        im = ax.pcolormesh(ds_ice.newlon, ds_ice.newlat, ice_year.isel(time = 0), vmin = 0, vmax = 100,
                           transform=ccrs.PlateCarree(), cmap=cmocean.cm.ice)
        ax.pcolormesh(ds_ice.newlon, ds_ice.newlat, test.isel(time = 0), vmin = 0, vmax = 1,
                           transform=ccrs.PlateCarree(), cmap = plt.cm.spring_r, alpha = 0.8)
        c = str(10*i+10) + "%"
        ax.text(0,-90, c, color='w', ha='center')
        if i == 7:
            new_ice = change_start_x(ds_ice, 0)
            new_ice_year = new_ice.siconc.isel(time = slice(188, 189))
            ax.contour(new_ice.newlon, new_ice.newlat, new_ice_year.isel(time = 0), levels=[90], colors='red', transform=ccrs.PlateCarree())
    
    fig.savefig(savename, dpi = 300)

def plot_sub_corr_plot(ax, typedf, df, col1, col2, col1_label, col2_label, sublabel, marker_dict, color_dict, labelfs = 10, loc = (0.93, 0.95)):
    dfa = pd.merge(typedf, df[['name', col1, col2]], on="name")
    dfa = dfa.dropna()
    plot_by_type(ax, dfa, marker_dict, color_dict)
    add_text_annotation(ax, dfa, t = 'number', expand_rate = 1, fs = 8)
    # ax1.set_xlim(-0.5, 10.5)
    modify_corr_ax(ax, col1_label, col2_label, fs = labelfs)
    ax.set_box_aspect(1)
    ax.tick_params(axis='both', labelsize=8)
    ax.text(loc[0], loc[1], sublabel, transform=ax.transAxes, fontsize=10, fontweight='bold')
    return ax
    
def plot_corr_4(typedf, df_pro, marker_dict, color_dict, loc=(-0.18, 0.95), savename = 'mean55.pdf'):
    fig = plt.figure(constrained_layout=True, figsize=(7, 7))

    axs = fig.subplots(nrows=2, ncols=2)
    ax = plot_sub_corr_plot(axs[0, 0], typedf, df_pro, 't_mean55', 'ice_mean55', 
                            "mean SST ($^{\circ}$C)", "mean SIC (%)", "a.",
                            marker_dict, color_dict, loc=loc)
    ax.legend(prop={'size': 9})
 
    ax = plot_sub_corr_plot(axs[0, 1], typedf, df_pro, 'hfds_m55', 'ice_mean55', 
                            "mean downward heatflux ($W/m^2$)", "mean SIC (%)",  "b.",
                            marker_dict, color_dict, loc=loc)
    
    ax = plot_sub_corr_plot(axs[1, 0], typedf, df_pro, 't_mean55', 'thick_m55', 
                            "mean SST ($^{\circ}$C)", "mean ice thickness (m)",  "c.",
                            marker_dict, color_dict, loc=loc)
 
    ax = plot_sub_corr_plot(axs[1, 1], typedf, df_pro, 'hfds_m55', 'mld_m55', 
                            "mean downward heatflux ($W/m^2$)", "mean MLD (m)",  "d.",
                            marker_dict, color_dict, loc=loc)
    plt.savefig(savename, bbox_inches='tight')


def plot_comparison(df_plot, df_plot2, marker_dict, color_dict, savename = 'comparison.pdf'):
    fig = plt.figure(constrained_layout=True, figsize=(7, 3.5))

    ax1, ax2 = fig.subplots(nrows=1, ncols=2)
    
    plot_by_type(ax1, df_plot, marker_dict, color_dict)
    df11 = df_plot[(df_plot['c_total'] >= 0.3) | (df_plot['p_total']  >= 0.3)]
    add_text_annotation(ax1, df11, t = 'number', fs = 8)
    ax1.plot([-1, 11],[-1, 11], color='black', linewidth=0.8, linestyle = "--", alpha = 0.3, zorder=1)
    ax1.set_xlim(-0.5, 10.5)
    ax1.set_ylim(-0.5, 10.5)
    modify_corr_ax(ax1, "Total convection region area ($10^6\ km^2$)", "Total polynya region area ($10^6\ km^2$)", equalbox = True, fs = 10)
    ax1.text(9.5, -0.1, "a.", fontsize = 12, fontweight = "bold")
    
    
    left_inset_ax = fig.add_axes([.30, .72, .25, .25], zorder = 1)
    left_inset_ax.set_xlim(-0.025, 0.25)
    left_inset_ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    df12 = df_plot[(df_plot['c_total'] < 0.3) & (df_plot['p_total']  < 0.3)]
    plot_by_type(left_inset_ax, df12, marker_dict, color_dict)
    add_text_annotation(left_inset_ax, df12, t = 'number', fs = 7, expand_rate = 2.5)
    modify_corr_ax(left_inset_ax, '', '', equalbox = True, fs = 8)
    
    
    plot_by_type(ax2, df_plot2, marker_dict, color_dict)
    add_text_annotation(ax2, df_plot2, t = 'number', fs = 8, expand_rate = 2.5)
    ax2.set_xlim(-5, 95)
    modify_corr_ax(ax2, "overlap ralative to convection (%)", "overlap ralative to polynya (%)", equalbox = True, fs = 10)
    
    handles, labels = ax1.get_legend_handles_labels()
    ax2.legend(handles, labels, prop={'size': 9})
    ax2.text(85, -0.1, "b.", fontsize = 12, fontweight = "bold")
    
    plt.savefig(savename, bbox_inches='tight')



def plot_corr_6(typedf, df_pro, marker_dict, color_dict, loc=(-0.22, 0.95), savename = 'corr_conv.pdf'):
    fig = plt.figure(constrained_layout=True, figsize=(9, 5.5))   
    axs = fig.subplots(nrows=2, ncols=3)
    
    area_lable = "mean convection area " +  r"($\times 10^6\ km^2$)"
    lfs = 8
    
    ax = plot_sub_corr_plot(axs[0, 0], typedf, df_pro,  't_mean_c', 'ice_mean_c',
                            "mean SST ($^{\circ}$C)","mean SIC (%)", "a.",
                            marker_dict, color_dict, labelfs = lfs, loc=loc)
    ax.legend(prop={'size': 9})
    
    ax = plot_sub_corr_plot(axs[0, 1], typedf, df_pro, 'hfds_mc', 'ice_mean_c', 
                            "mean downward heatflux ($W/m^2$)","mean SIC (%)", "b.",
                            marker_dict, color_dict, labelfs = lfs, loc=loc)
    ax = plot_sub_corr_plot(axs[0, 2], typedf, df_pro, 'ice_mean_c', 'c_mean', 
                            "mean SIC (%)", area_lable,  "c.",
                            marker_dict, color_dict, labelfs = lfs, loc=loc)
    
    ax = plot_sub_corr_plot(axs[1, 0], typedf, df_pro, 'tann_c', 'c_mean', 
                            "annual mean SST ($^{\circ}$C)", area_lable , "d.",
                            marker_dict, color_dict, labelfs = lfs, loc=loc)
    
    ax = plot_sub_corr_plot(axs[1, 1], typedf, df_pro, 'mld_mc', 'c_mean', 
                            "mean september MLD (m)", area_lable,  "e.",
                            marker_dict, color_dict, labelfs = lfs, loc=loc)
    ax = plot_sub_corr_plot(axs[1, 2], typedf, df_pro, 'hfds_mc', 'c_mean', 
                            "mean downward heatflux ($W/m^2$)", area_lable,  "f.",
                            marker_dict, color_dict, labelfs = lfs, loc=loc)

    plt.savefig(savename, bbox_inches='tight')

def add_subplot_num(ax, t, y, x = -1.4):
    ax.text(x, y, t, fontweight="bold")

def add_bars(ax, i, df, color_dict, phn = 9):
    df_plot_i = df[['name', 'type', 'number', df.columns.tolist()[i], df.columns.tolist()[i+1]]]
    ph = np.max(np.abs(df_plot_i[df_plot_i.columns.tolist()[3]]))/phn
    for index, row in df_plot_i.iterrows():
        ax.bar(
            str(row['number']), row[df.columns.tolist()[i]],
            edgecolor = 'k', color=color_dict[row['type']], zorder=3
        )
        if row[df.columns.tolist()[i+1]] < 0.005:
            if row[df.columns.tolist()[i]] > 0:
                ht = row[df.columns.tolist()[i]] + ph
            else:
                ht = row[df.columns.tolist()[i]] - ph
            ax.scatter(str(row['number']), ht, marker = "x", color = 'k', s = 10)
    # ax.spines[:].set_visible(False)

def add_bar_plot(ax0, df_plot4, color_dict, i, ylabel, pos, yls, yticks, dy, lcoords, phn=9):
    twin = ax0.twinx()
    add_bars(twin, i, df_plot4, color_dict, phn=phn)
    twin.spines[:].set_visible(False)
    twin.spines.left.set_visible(True)
    # twin.spines.left.set_color(p2.get_color())
    # twin.spines.left.set_position(("axes", -0.085))
    if pos == 'left':
        twin.tick_params(left=True, labelleft=True, right=False, labelright=False, direction='inout')
    elif pos == 'right':
        twin.tick_params(left=False, labelleft=False, right=True, labelright=True, direction='inout')
    # twin.tick_params(axis='both', which='major', labelsize=6)
    # twin.tick_params(axis='y', colors=p2.get_color(), labelrotation=90)
    twin.set_ylabel(ylabel)
    # twin.yaxis.set_label_coords(-0.15, 0.5)
    twin.yaxis.set_label_position(pos)
    twinyticks = twin.get_yticklabels()
    twin.set_ylim(yls[0], yls[1])
    twin.spines.left.set_bounds(yticks[0], yticks[1])
    twin.set_yticks(np.arange(yticks[0], yticks[1]+dy, dy))
    twin.yaxis.set_label_coords(lcoords[0], lcoords[1])
    twin.axhline(0, color = 'k', linewidth = 0.5)

def plot_deep_shallow(df_plot4, color_dict, savename = 'deepshallow.pdf'):
    fig = plt.figure(constrained_layout=True, figsize=(7, 7))

    filtered_colors = {k: color_dict[k] for k in set(df_plot4['type'])}
    add_type_color_legend(fig, filtered_colors, '', bta = (0.85, 0.99), n_cols = 6, fs=8)
    
    ax0 = fig.subplots()
    ax0.bar([str(num) for num in df_plot4['number'].tolist()], [0 for num in df_plot4['number'].tolist()])
    ax0.set_yticks(np.arange(0, 150, 10))
    ax0.set_ylim(0, 150)
    ax0.tick_params(axis='both', direction='in')
    ax0.tick_params(axis='y', which='both', left=False, labelleft=False)
    ax0.grid(True)
    add_subplot_num(ax0, 'a.', 136)
    add_subplot_num(ax0, 'b.', 116)
    add_subplot_num(ax0, 'c.', 106)
    add_subplot_num(ax0, 'd.', 88)
    add_subplot_num(ax0, 'e.', 73)
    add_subplot_num(ax0, 'f.', 49)
    add_subplot_num(ax0, 'g.', 27)
    add_subplot_num(ax0, 'h.', 13)
    
    add_bar_plot(ax0, df_plot4, color_dict, 13, "MLD (m)", "left", (0, 60000), (0, 8000), 2000, (-0.08, 0.08), phn=13)
    add_bar_plot(ax0, df_plot4, color_dict, 11, "downward heatflux\n(W/m$^2$)", "right", (-1400, 1600), (-200, 100), 50, (1.09, 0.44), phn=7)
    add_bar_plot(ax0, df_plot4, color_dict, 9, "ice thickness (m)", "left", (-1.6, 4.4), (-0.2, 0.6), 0.2, (-0.08, 0.29), phn=11)
    add_bar_plot(ax0, df_plot4, color_dict, 7, "SIC (%)", "right", (-90, 360), (-30, 10), 10, (1.08, 0.18), phn=7)
    add_bar_plot(ax0, df_plot4, color_dict, 3, "SST ($^{\circ}$)", "left", (-16, 14), (-0.5, 2), 0.5, (-0.08, 0.55), phn=4)
    add_bar_plot(ax0, df_plot4, color_dict, 15, "annual SST ($^{\circ}$)", "right", (-20, 10), (-1, 1.5), 0.5, (1.08, 0.670), phn=3)
    add_bar_plot(ax0, df_plot4, color_dict, 5, "SSS", "left", (-6.6, 2.4), (-0, 0.8), 0.2, (-0.08, 0.80), phn=5)
    add_bar_plot(ax0, df_plot4, color_dict, 17, "annual SSS", "right", (-7.8, 1.2), (0, 0.8), 0.2, (1.07, 0.92), phn=5)
    
    plt.savefig(savename, bbox_inches='tight')




def main():
    # filter some warning messages
    import warnings
    warnings.filterwarnings("ignore")
    
    from mycolors import PALETTE as color_palette
    from mycolors import MARKERS as marker_list

    p_sos = '../../SO_data/data_sos/'
    p_sos_ann = '../../SO_data/data_sos_ann/'
    p_sst = '../../SO_data/data_sst/'
    p_sst_ann = '../../SO_data/data_sst_ann/'
    p_hfds = '../../SO_data/data_hfds/'
    p_ice = '../../SO_data/data_siconc_w_area/'
    p_mld = '../../SO_data/data_mld/'
    p_mlotst = '../../SO_data/data_mlotst/'
    p_thick = '../../SO_data/data_thick/'
    p_polynya = '../../SO_data/data_polynya_mean/'

    name = 'GFDL-CM4'
    data_siconc = openpickle(name, p_ice)
    data_sithick = openpickle(name, p_thick)
    data_mld = openpickle(name, p_mld)
    plot_methods(data_mld, data_siconc, data_sithick)
    
    ds_ice_id = openpickle("NorESM2-MM", p_ice)
    plot_detection(ds_ice_id)
    
    df_pro = pd.read_csv('properties.csv')
    dfs = df_pro[['name', 'c_total', 'p_total']]
    
    df_cross = pd.read_csv('CrossRate.csv')
    dfs2 = df_cross[['name', 'c_cross_rate', 'p_cross_rate']]
    dfs3 = df_cross[['name', 'p_events', 'c_events', 'p_cross_c_events']]
    
    datapd = pd.read_csv('List_model.csv')
    typedf = datapd[['source_id', 'type']]
    typedf['number'] = datapd.index + 1
    typedf = typedf.rename(columns={'source_id': 'name'})
    
    color_dict = create_new_color_dict(color_palette, typedf['type'].unique())
    marker_dict = create_new_marker_dict(marker_list, typedf['type'].unique())
    
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

    df_dpsp = pd.read_csv('properties_pcdiff_ttest.csv')
    df_plot4 = pd.merge(typedf, df_dpsp, on="name")
    df_plot4 = df_plot4[(df_plot4['ice_pcnc'].notna())]

    plot_corr_4(typedf, df_pro, marker_dict, color_dict)
    plot_corr_6(typedf, df_pro, marker_dict, color_dict)
    plot_comparison(df_plot, df_plot2, marker_dict, color_dict)
    plot_deep_shallow(df_plot4, color_dict)


if __name__ == "__main__":
    main()
