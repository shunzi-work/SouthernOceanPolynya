###############
# 
# # Python script for plotting polynya
# 
###############

from myfunctions import *
import gc
import matplotlib.pyplot as plt
from IcePlot import create_new_color_dict, add_subplot_icepolynya, add_cbars, add_type_color_legend

def plot_polynya_maps(datapd, p_ice, p_polynya, colors, typename='type_ice', titleonly=True, resotext=False, lighten =1, figsize=(6.5,7.5)):
    fig = plt.figure(figsize=figsize)
    n = 0
    color_dict = create_new_color_dict(colors, datapd, typename, lighten = lighten)
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        t = datapd.at[i, typename]
        dsice = openpickle(name, p_ice)
        plt_icemax = dsice.siconc.max('time')
        ds_polynya = openpickle(name, p_polynya)
        timelength = len(ds_polynya.time)
        plt_polynya = ds_polynya.count('time').where(ds_polynya.count('time')>0)/timelength
        
        pltx = dsice.newlon
        plty = dsice.newlat

        n = n + 1
        if resotext:
            reso_text = str(datapd.at[i, 'resolution'])
        else:
            reso_text = None
        im, im2 = add_subplot_icepolynya(fig, n, name, color_dict[t], titleonly, reso_text, pltx,plty,plt_icemax,plt_polynya, timelength)
        gc.collect()
    add_cbars(fig, im, im2)
    add_type_color_legend(fig, list(color_dict.keys()),  list(color_dict.values()), 'sea ice module types')
    figsavename = 'Figures/polynya_maps.png'
    fig.savefig(figsavename, dpi = 300)




def plot_convection_maps(datapd, p_mld, colors, typename='type', titleonly=True, resotext=True, lighten=1, figsize=(6.5,7.5)):
    fig = plt.figure(figsize=figsize)
    n = 0
    color_dict = create_new_color_dict(colors, datapd, typename, lighten = lighten)
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        t = datapd.at[i, typename]
        dsice = openpickle(name, p_ice)
        plt_icemax = dsice.siconc.max('time')
        ds_polynya = openpickle(name, p_polynya)
        timelength = len(ds_polynya.time)
        plt_polynya = ds_polynya.count('time').where(ds_polynya.count('time')>0)/timelength
        
        pltx = dsice.newlon
        plty = dsice.newlat

        n = n + 1
        if resotext:
            reso_text = str(datapd.at[i, 'resolution'])
        else:
            reso_text = None
        im, im2 = add_subplot_icepolynya(fig, n, name, color_dict[t], titleonly, reso_text, pltx,plty,plt_icemax,plt_polynya, timelength)
        gc.collect()
    add_cbars(fig, im, im2)
    add_type_color_legend(fig, list(color_dict.keys()),  list(color_dict.values()), 'sea ice module types')
    figsavename = 'Figures/polynya_maps.png'
    fig.savefig(figsavename, dpi = 300)