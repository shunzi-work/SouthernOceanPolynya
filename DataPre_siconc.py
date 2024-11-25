###############
# 
# # Python script for CMIP6 siconc data preprocessing 
# 
###############

from myfunctions import *

def get_new_siconc_dataset(datapd, p_save, p_nc, selected_month, southlat, newx=False):
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        print("{} {}".format(i, name), end = ' ')
        
        # check if the resulting file exists
        if ispickleexists(name, p_save):
            print("    [o] data exist.")
            continue
    
        # check if siconc data is avalaible online
        if pd.isna(datapd.at[i, 'zstore_siconc']):
            ds = read_siconc(p_nc, name, selected_month)
        else:
            ds = open_from_cloud(datapd.at[i, 'zstore_siconc'])
            ds = select_month(ds, selected_month)
    
        icedata = ds.siconc
        icedata = icedata.where(icedata>=0).where(icedata<=100)
    
        nolatlon = False
        newlatlon = False
        
        if name == 'NESM3':
            newlatlon = ('lat', 'lon')
            
        dlat, dlon = get_latlon(datapd, i, icedata.isel(time=0), newlatlon, nolatlon)
    
        if dlat.values.flatten()[~np.isnan(dlat.values.flatten())][0]>0:
            icedata = flip_y(icedata)
            dlat, dlon = get_latlon(datapd, i, icedata.isel(time=0), newlatlon, nolatlon)
        
        if ('area' in ds) or ('areacello' in ds):
            if 'area' in ds:
                if 'areacello' not in ds:
                    ds = ds.rename({'area':'areacello'})
                    
            dsg_data = ds.areacello
            if 'time' in dsg_data.dims:
                dsg_data = dsg_data.isel(time = 0)
            dsg_data = dsg_data.where(dsg_data>0).where(dsg_data<1e30)
            dlat_g, dlon_g = dlat, dlon
        else:
            if name == 'GISS-E2-2-H':
                nameg = 'GISS-E2-1-H'
            elif name == 'UKESM1-1-LL':
                nameg = 'UKESM1-0-LL'
            else:
                nameg = name
    
            if name == 'NESM3':
                varname = 'areacelli'
            else:
                varname = 'areacello'
    
            dsg = read_areacello(p_nc, nameg, varname) 
            dsg_data = dsg[varname]
            dsg_data = dsg_data.where(dsg_data>0).where(dsg_data<1e30)
            
            if name == 'CAS-ESM2-0': # CAS-ESM2-0 cell area grid has different name
                nolatlon = True
                newlatlon = ('lat', 'lon')
            
            dlat_g, dlon_g = get_latlon(datapd, i, dsg_data, newlatlon, nolatlon)
    
            if dlat_g.values.flatten()[~np.isnan(dlat_g.values.flatten())][0]>0:
                dsg_data = flip_y(dsg_data)
                dlat_g, dlon_g = get_latlon(datapd, i, dsg_data, newlatlon, nolatlon)
    
        if dlat.shape != dlat_g.shape: # if not the same shape
            if name in ['CMCC-CM2-SR5', 'CMCC-ESM2']:
                dsg_data = dsg_data.isel({dsg_data.dims[-1] : slice(1, len(dsg_data[dsg_data.dims[-1]])-1)})
                dsg_data = dsg_data.isel({dsg_data.dims[-2] : slice(0, len(dsg_data[dsg_data.dims[-2]])-1)})
                dsg_data = copy_xy(icedata, dsg_data)
                dlat_g, dlon_g = get_latlon(datapd, i, dsg_data, newlatlon, nolatlon)
            elif name in ['NorESM2-MM', 'NorESM2-LM']:
                dsg_data = dsg_data.isel({dsg_data.dims[-2] : slice(0, len(dsg_data[dsg_data.dims[-2]])-1)})
                dlat_g, dlon_g = get_latlon(datapd, i, dsg_data, newlatlon, nolatlon)
            else:
                print('    Skip.')
                continue
                
        if dlat.shape == dlat_g.shape: # same shape
            if (np.nanmax(np.abs(dlat.values - dlat_g.values)) < 10e-4) and (np.nanmax(np.abs(dlon.values - dlon_g.values)) < 10e-4):
                if np.isnan(dlat).any() or np.isnan(dlon).any(): # if coords in ice data is not complete
                    if np.isnan(dlat_g).any() or np.isnan(dlon_g).any(): # if coords in area data is not complete
                        newdlon, newdlat = newxy_fmissingxy(dlon, dlat)
                        ds_siconc = drop_coords(icedata)
                        new_ds = create_new_ds(ds_siconc, dsg_data, newdlat, newdlon)
                        print('    No complete coords, fill missing coords ...', end = ' ')
                    else:
                        ds_siconc = copy_xy(dsg_data, icedata)
                        ds_siconc = drop_coords(ds_siconc)
                        new_ds = create_new_ds(ds_siconc, dsg_data, dlat_g, dlon_g)
                        print('    No complete coords, use coords from area data ...', end = ' ')
                else:
                    new_ds = create_new_ds(icedata, dsg_data, dlat, dlon)
                    print('    Coords match ...', end = ' ')
            else:
                if name in ['CESM2-WACCM-FV2', 'NorESM2-MM', 'NorESM2-LM','FGOALS-g3', 'CAS-ESM2-0']:
                    ds_siconc = copy_xy(dlat_g, icedata)
                    ds_siconc = drop_coords(ds_siconc)
                    new_ds = create_new_ds(ds_siconc, dsg_data, dlat_g, dlon_g)
                    print('    Not match, use coords from area data ...', end = ' ')
                else:
                    if pd.isna(datapd.at[i, 'latname']):
                        dsg_data = calculate_area_xy(ds)
                        new_ds = create_new_ds(icedata, dsg_data, dlat, dlon)
                        print('    No area data, calculate area ...', end = ' ')                   
                    else:
                        print('    Skip.')
                        continue
        
        slat = new_ds['newlat'].where(new_ds['newlat']<=southlat, drop = True)
        new_ds_south = new_ds.sel({slat.dims[0] : slat[slat.dims[0]]})
        
        if new_ds_south.siconc.isel(time = 0)[0,0] == 0:
            ### ['GISS-E2-1-H', 'GISS-E2-2-H', 'INM-CM4-8']
            new_ds_south = set_land_to_nan(new_ds_south)      

        ## get rid of blank columns     
        while np.isnan(new_ds_south.siconc.isel({new_ds_south.siconc.dims[-1]:0}).isel(time=0)).all():
            new_ds_south = new_ds_south.isel({new_ds_south.siconc.dims[-1]:slice(1, None)})
        while np.isnan(new_ds_south.siconc.isel({new_ds_south.siconc.dims[-1]:-1}).isel(time=0)).all():
            new_ds_south = new_ds_south.isel({new_ds_south.siconc.dims[-1]:slice(0, -1)})
            
        if np.isnan(new_ds_south.siconc.isel(time = 0)[-1,-1]):
            ### ['E3SM-2-0', 'E3SM-2-0-NARRM']
            new_ds_south = set_ocean_to_zero(new_ds_south)
            
        new_ds_south = new_ds_south.where(new_ds_south['newlat'] <= southlat, drop = True)
        if newx: 
            new_ds_south = change_start_x(new_ds_south, newx)   
        new_ds_south.load()
        savepickle(name, p_save, new_ds_south)
        print('[*] Saved.')

def count_area_diff_thresholds(datapd, ice_thresholds, area_threshold, p_save, datap0, re=False):
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        print("{} {} ....".format(i, name), end = ' ')
    
        if ispickleexists(name, p_save):
            print("[o] exist.")
            continue
            
        if ispickleexists(name, datap0):
            ds = openpickle(name, datap0)
            area_count = []
            for ice_threshold in ice_thresholds:
                area_count.append(count_polynya_area(ds, ice_threshold, area_threshold, re=re))
            savepickle(name, p_save, area_count)
            print("[*] saved.")

def main():
    # filter some warning messages
    import warnings
    warnings.filterwarnings("ignore")
    
    datapd = pd.read_csv('List57.csv')
    p_save = '../../SO_data/data_siconc_w_area/'
    p_nc = '../../../data/model/CMIP6/'
    selected_month = 9
    southlat = -40 
    newx = 135

    print('Start siconc data preprocessing ...')
    get_new_siconc_dataset(datapd, p_save, p_nc, selected_month, southlat, newx)

    print('Finish siconc data preprocessing.')

    print()
    print('Start polynya detecting ...')
    print('Count polynya area using different thresholds ...')

    ice_thresholds = np.arange(0, 100, step=1)
    area_threshold = [10000, 1000000]
    # count_save = '../../SO_data/data_polynya_count/'
    count_save = '../../SO_data/data_polynya_count_re/'
    count_area_diff_thresholds(datapd, ice_thresholds, area_threshold, count_save, p_save, re=2)

if __name__ == "__main__":
    main()