###############
# 
# # Python script for plotting polynya
# 
###############

from myfunctions import *
from PlotFunctions import *
from mycolors import PALETTE as color_palette


import gc

def main():
    # filter some warning messages
    import warnings
    warnings.filterwarnings("ignore")
    
    datapd = pd.read_csv('List_model.csv')
    p_ice = '../../SO_data/data_siconc_w_area/'
    p_count = '../../SO_data/data_polynya_count/'
    p_polynya = '../../SO_data/data_polynya_mean/'
    p_mld = '../../SO_data/data_mld/'
    p_mlotst = '../../SO_data/data_mlotst/'

    # print('Plot histogram ...')
    # plot_ice_hist(datapd, p_ice, p_count, 'PolynyaCounts.pdf')
    # gc.collect()

    # print('Plot polynya maps (different thresholds) ...')
    # ice_thresholds = np.arange(0, 100, step=1)
    # savename_pre = 'Figures/polynya_maps_'
    # for ice_threshold in ice_thresholds:
    #     sn = savename_pre + str(ice_threshold) + '.png'
    #     if os.path.exists(sn):
    #         continue
    #     plot_polynya_maps(datapd, p_ice, color_palette, ice_threshold, savename_pre = savename_pre)
    #     gc.collect()
    
    print('Plot polynya maps (from polynya file) ...')
    plot_polynya_maps_from_polynya_data(datapd, p_ice, p_polynya, color_palette)
    gc.collect()

    print('Plot convection maps ...')
    plot_convection_maps(datapd, p_mld, p_mlotst, color_palette)
    gc.collect()


if __name__ == "__main__":
    main()
