###############
# 
# # Python script for CMIP6 siconc data preprocessing 
# 
###############

from myfunctions import *
import gc


def save_new_dataset(datapd, p_save, p_nc, selected_month, southlat, dataname, newx=False):
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        print("{} {}".format(i, name), end = '...')
        if ispickleexists(name, p_save):
            print("[o] data exist.")
            continue
        new_ds_south = get_new_dataset(datapd.iloc[i], p_nc, selected_month, southlat, dataname, newx=newx)
        if isinstance(new_ds_south, xr.Dataset):
            new_ds_south = new_ds_south.load()
            savepickle(name, p_save, new_ds_south)
            print("[*] Saved.")
            gc.collect()
        tempfile_path = p_nc + '*' + '_temp.nc'
        matching_files_temp = glob.glob(tempfile_path)
        if len(matching_files_temp)>0:
            for mf in matching_files_temp:
                os.remove(mf)

def get_new_dataset(data_info, p_nc, selected_month, southlat, dataname, newx=False):
    # check if mlotst data is avalaible on cloud
    datastorename = 'zstore_' + dataname 
    name = data_info['source_id']
    if pd.isna(data_info[datastorename]):
        try:
            ds = read_nc_files(p_nc, name, selected_month, dataname)
        except Exception as e:
            print(e)
            return None
    else:
        if data_info[datastorename][0:2] == 'gs':
            ds = open_from_cloud(data_info[datastorename]) 
            if selected_month>0:
                ds = select_month(ds, selected_month)
        else:
            datafile_paths = data_info[datastorename] + '*'
            matching_files = glob.glob(datafile_paths)
            if selected_month>0:
                ds = open_nc_month(matching_files, selected_month)
            else:
                ds = open_nc(matching_files)
    ds = ds.sortby('time')

    dataarray = ds[dataname]
    if dataname == 'siconc':
        dataarray = dataarray.where(dataarray>=0).where(dataarray<=100)
    elif dataname in ['so', 'sos']:
        dataarray = dataarray.where(dataarray>=10).where(dataarray<50)
        # E3SM missing value is 1
    elif dataname in ['thetao', 'tos']:
        dataarray = dataarray.where(dataarray!=0).where(dataarray<100)
    elif dataname == 'hfds':
        if name == 'CAS-ESM2-0':
            # https://errata.ipsl.fr/static/view.html?uid=34c23dd3-2aac-f635-8538-eaf110e83611
            # All the hfds data are saved as K/s instead of W/m2, they 
            # should be multiplied by factor of 4.1*10^7 for unit 
            # conversion from K/s to W/m2
            dataarray = dataarray * 4.1e7
        if name == "SAM0-UNICON":
            # No errata info, but hfds sign should be flipped
            dataarray = dataarray * (-1)

    nolatlon = False
    newlatlon = False

    if (dataname in ['siconc', 'sithick']) and (name == 'NESM3'):
        newlatlon = ('lat', 'lon')

    if name == 'CAS-ESM2-0':
        if dataname in ['mlotst', 'thetao', 'so', 'hfds', 'tos', 'sos']:
            # CAS-ESM2-0 sea water dataset grid different from sea ice dataset grid
            # sea water: lat, lon;  sea ice: i, j, longitude, latitude
            nolatlon = True
            newlatlon = ('lat', 'lon')
        else:
            ## 'CAS-ESM2-0' sea water lon 0~359 ; ice lon 1~360 
            ## shift ice data -> 0 ~ 359
            dataarray = shift_x(dataarray)

    dlat, dlon = get_latlon(data_info, dataarray.isel(time=0), newlatlon, nolatlon)
    if dlat.values.flatten()[~np.isnan(dlat.values.flatten())][0]>0:
        dataarray = flip_y(dataarray)
        dlat, dlon = get_latlon(data_info, dataarray.isel(time=0), newlatlon, nolatlon)
        
    if ('area' in ds) or ('areacello' in ds):
        print('areadata exist...', end = '')
        if 'area' in ds:
            if 'areacello' not in ds:
                ds = ds.rename({'area':'areacello'})
                print('rename area to areacello...', end = '')
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
        if (name == 'NESM3') and (dataname in ['siconc', 'sithick']):
            varname = 'areacelli'
        else:
            varname = 'areacello'
        dsg = read_areacello(p_nc, nameg, varname, data_info) 
        dsg_data = dsg[varname]
        dsg_data = dsg_data.where(dsg_data>0).where(dsg_data<1e20)
        
        if name == 'CAS-ESM2-0': # CAS-ESM2-0 cell area grid has different name
            nolatlon = True
            newlatlon = ('lat', 'lon')

        if name == 'NESM3':
            newlatlon = ('lat', 'lon')
        
        dlat_g, dlon_g = get_latlon(data_info, dsg_data, newlatlon, nolatlon)

        if dlat_g.values.flatten()[~np.isnan(dlat_g.values.flatten())][0]>0:
            dsg_data = flip_y(dsg_data)
            dlat_g, dlon_g = get_latlon(data_info, dsg_data, newlatlon, nolatlon)
     
    if dlat.shape != dlat_g.shape: # if not the same 
        print('!!!! not the same', end = '')
        if dataname in ['siconc', 'sithick']:
            if name in ['CMCC-CM2-SR5', 'CMCC-ESM2']:
                dsg_data = dsg_data.isel({dsg_data.dims[-1] : slice(1, len(dsg_data[dsg_data.dims[-1]])-1)})
                dsg_data = dsg_data.isel({dsg_data.dims[-2] : slice(0, len(dsg_data[dsg_data.dims[-2]])-1)})
                dsg_data = copy_xy(dataarray, dsg_data)
                dlat_g, dlon_g = get_latlon(data_info, dsg_data, newlatlon, nolatlon)
            elif name in ['NorESM2-MM', 'NorESM2-LM']:
                dsg_data = dsg_data.isel({dsg_data.dims[-2] : slice(0, len(dsg_data[dsg_data.dims[-2]])-1)})
                dlat_g, dlon_g = get_latlon(data_info, dsg_data, newlatlon, nolatlon)
            else:
                print('......Skip.')
                return None
        # elif dataname == 'hfds':

    if (dlat.shape == dlat_g.shape) and (name not in ['KIOST-ESM']): # same shape
        if (np.nanmax(np.abs(dlat.values - dlat_g.values)) < 10e-4) and (np.nanmax(np.abs(dlon.values - dlon_g.values)) < 10e-4):
            if np.isnan(dlat).any() or np.isnan(dlon).any(): # if coords in ice data is not complete
                if np.isnan(dlat_g).any() or np.isnan(dlon_g).any(): # if coords in area data is not complete
                    newdlon, newdlat = newxy_fmissingxy(dlon, dlat)
                    dataset_final = drop_coords(dataarray)
                    new_ds = create_new_ds(dataset_final, dsg_data, newdlat, newdlon, dataname)
                    print('    No complete coords, fill missing coords ...', end = ' ')
                else:
                    dataset_final = copy_xy(dsg_data, dataarray)
                    dataset_final = drop_coords(dataset_final)
                    new_ds = create_new_ds(dataset_final, dsg_data, dlat_g, dlon_g, dataname)
                    print('    No complete coords, use coords from area data ...', end = ' ')
            else:
                dataset_final = copy_xy(dsg_data, dataarray)
                dataset_final = drop_coords(dataset_final)
                new_ds = create_new_ds(dataset_final, dsg_data, dlat, dlon, dataname)
                print('    Coords match ...', end = ' ')
        else:
            if name in ['CESM2-WACCM-FV2', 'NorESM2-MM', 'NorESM2-LM','FGOALS-g3', 'CAS-ESM2-0']:
                dataset_final = copy_xy(dsg_data, dataarray)
                dataset_final = drop_coords(dataset_final)
                new_ds = create_new_ds(dataset_final, dsg_data, dlat_g, dlon_g, dataname)
                print('    Not match, use coords from area data ...', end = ' ')
            else:
                if pd.isna(data_info['latname']):
                    dsg_data = calculate_area_xy(ds, dataname)
                    new_ds = create_new_ds(dataarray, dsg_data, dlat, dlon, dataname)
                    print('    No area data, calculate area ...', end = ' ')                   
                else:
                    dsg_data = calculate_area_latlon(ds, data_info)
                    new_ds = create_new_ds(dataarray, dsg_data, dlat, dlon, dataname)
                    print('    No area data, calculate area (latlon) ...', end = ' ')     
    else:
        if name in ['KIOST-ESM']:
            dsg_data = calculate_area_xy(ds, dataname)
            new_ds = create_new_ds(dataarray, dsg_data, dlat, dlon, dataname)
            print('    No area data, calculate area ...', end = ' ')                   
        elif name in ['GISS-E2-1-H', 'GISS-E2-2-H']:
            dataarray = regrid_based_on_dsgxy(dataarray, dsg_data, data_info)
            new_ds = create_new_ds(dataarray, dsg_data, dlat_g, dlon_g, dataname)
            print('    regriding ...', end = ' ')   
        else:
            print("HERE?")
            return None

    # south (first use select method to avoid nans)
    slat = new_ds['newlat'].where(new_ds['newlat']<=southlat, drop = True)
    new_ds_south = new_ds.sel({slat.dims[0] : slat[slat.dims[0]]})
    
    if (dataname == 'siconc') and (new_ds_south[dataname].isel(time = 0)[0,0] == 0):
        ### ['GISS-E2-1-H', 'GISS-E2-2-H', 'INM-CM4-8']
        new_ds_south = set_land_to_nan(new_ds_south) 

    ## get rid of blank columns    
    # if dataname in ['siconc','mlotst','hfds', 'so', 'sos', 'tos', 'thetao']: 
    while np.isnan(new_ds_south[dataname].isel({new_ds_south[dataname].dims[-1]:0}).isel(time=0)).all():
        new_ds_south = new_ds_south.isel({new_ds_south[dataname].dims[-1]:slice(1, None)})
    while np.isnan(new_ds_south[dataname].isel({new_ds_south[dataname].dims[-1]:-1}).isel(time=0)).all():
        new_ds_south = new_ds_south.isel({new_ds_south[dataname].dims[-1]:slice(0, -1)})
    
    if (dataname == 'siconc') and (np.isnan(new_ds_south[dataname].isel(time = 0)[-1,-1])):
        ### ['E3SM-2-0', 'E3SM-2-0-NARRM']
        new_ds_south = set_ocean_to_zero(new_ds_south)

    southnewlat = new_ds_south.newlat.copy()
    southnewlon = new_ds_south.newlon.copy()
        
    new_ds_south = new_ds_south.where(new_ds_south['newlat'] <= southlat, drop = True)
    new_ds_south['newlat'] = southnewlat
    new_ds_south['newlon'] = southnewlon
    if newx: 
        new_ds_south = change_start_x(new_ds_south, newx)

    if selected_month!=0:
        if (name in ['CAS-ESM2-0']) and (dataname in ['siconc', 'sivol', 'sithick']):
            #CAS-ESM2-0 with one year less in mld data than in ice data
            new_ds_south = new_ds_south.isel(time = slice(0, -1))                
        if len(new_ds_south.time) > 500:
            # if (name in ['SAM0-UNICON']) and (dataname == 'hfds'):
                # hfds dataset on NCAR has only 699 year (although with 700 years in names), no year 290
                # new_ds_south = new_ds_south
            if (name in ['ACCESS-ESM1-5']) and (dataname in ['tos', 'sos']):
                new_ds_south = new_ds_south.isel(time = slice(-600, -100))
            else:
                new_ds_south = new_ds_south.isel(time = slice(-500, None))
    return new_ds_south

def count_area_diff_thresholds(datapd, ice_thresholds, area_threshold, buffering, p_save, datap0, re=False):
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        print("{} {} ....".format(i, name), end = ' ')
    
        if ispickleexists(name, p_save):
            print("[o] exist.")
            continue
            
        if ispickleexists(name, datap0):
            ds = openpickle(name, datap0)
            flood_points = [(0,0)]
            if name == "MRI-ESM2-0":
                flood_points = [(0,0), (0,40), (0,100), (0,200)]
            area_count = []
            for ice_threshold in ice_thresholds:
                a0 = count_polynya_area(ds, ice_threshold, area_threshold, flood_points, buffering, re=re)
                area_count.append(a0)
            savepickle(name, p_save, area_count)
            print("[*] saved.")
            gc.collect()

def polynya_detecting_mean(datapd, p_save, p_ice, area_threshold):
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        print("{} {} ....".format(i, name), end = ' ')
                
        dssiconc = openpickle(name, p_ice)
        # ice_mean = (dssiconc.siconc*dssiconc.areacello).sum()/dssiconc.areacello.where(dssiconc.siconc>=0).sum()            
        ice_mean_not0 = (dssiconc.siconc.where(dssiconc.siconc.max("time")>0)*dssiconc.areacello).sum()/dssiconc.areacello.where(dssiconc.siconc.max("time")>0).sum()/len(dssiconc.time)

        flood_points = [(0,0)]
        if name == "MRI-ESM2-0":
            flood_points = [(0,0), (0,40), (0,100), (0,200)]

        if not ispickleexists(name, p_save):
            mask = detect_polynya(dssiconc.siconc, dssiconc.areacello, ice_mean_not0.values.item(), area_threshold, flood_points, buffering = 0)
            savepickle(name, p_save, mask)
            
        print("[*] saved.")
        gc.collect()


def polynya_detecting_fixed(datapd, p_save, p_ice, area_threshold, fixed_num):
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']

        print("{} {} ... ".format(i, name), end = '')
        if ispickleexists(name, p_save):
            print("[o] exist.")
            continue
        
        dssiconc = openpickle(name, p_ice)

        flood_points = [(0,0)]
        if name == "MRI-ESM2-0":
            flood_points = [(0,0), (0,40), (0,100), (0,200)]
        mask1 = detect_polynya(dssiconc.siconc, dssiconc.areacello, fixed_num, area_threshold, flood_points, buffering = 0)
        savepickle(name, p_save, mask1)

        print("[*] saved.")
        gc.collect()

def save_new_dataset_sithickness(datapd, p_save, p_nc, selected_month, southlat, newx=False):
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        print("{} {}".format(i, name), end = '...')
        if ispickleexists(name, p_save):
            print("[o] data exist.")
            continue
        new_ds_south = get_new_dataset(datapd.iloc[i], p_nc, selected_month, southlat, 'sithick', newx=newx)
        if not isinstance(new_ds_south, xr.Dataset):
            new_ds_south = get_new_dataset(datapd.iloc[i], p_nc, selected_month, southlat, 'sivol', newx=newx)
        if isinstance(new_ds_south, xr.Dataset):
            # if (name in ['SAM0-UNICON', 'CAS-ESM2-0']) and (dataname == 'siconc'):
            #     # SAM0-UNICON, CAS-ESM2-0 with one year less in mld data than in ice data
            #     new_ds_south = new_ds_south.isel(time = slice(0, -1))
            new_ds_south = new_ds_south.load()
            savepickle(name, p_save, new_ds_south)
            print("[*] Saved.")
            gc.collect()

def main():
    # filter some warning messages
    import warnings
    warnings.filterwarnings("ignore")
    
    datapd = pd.read_csv('List_model.csv')
    p_ice = '../../SO_data/data_siconc_w_area/'
    p_sith = '../../SO_data/data_thick/'

    p_nc = '../../data/CMIP6/'
    selected_month = 9
    southlat = -40 
    newx = 135
    dataname = 'siconc'

    print('Start siconc data preprocessing ...')
    # get_new_siconc_dataset(datapd, p_ice, p_nc, selected_month, southlat, newx)
    save_new_dataset(datapd, p_ice, p_nc, selected_month, southlat, dataname, newx=newx)
    print('Finish siconc data preprocessing.')

    print()
    print('Start polynya detecting ...')
    print('Count polynya area using different thresholds ...')
    ice_thresholds = np.arange(0, 100, step=1)
    area_threshold = [0, 2000]
    buffering = 15
    count_save = '../../SO_data/data_polynya_count/'
    count_area_diff_thresholds(datapd, ice_thresholds, area_threshold, buffering, count_save, p_ice)
    print('Finish polynya counting.')

    print()
    print('Start polynya detecting based on the ice threshold (mean SIC)')
    p_polynya_save = '../../SO_data/data_polynya_mean/'
    polynya_detecting_mean(datapd, p_polynya_save, p_ice, area_threshold)

    # print()
    # print("Start polynya detecting based on the a fixed thershold...")
    # p_polynya_save_fix = '../../SO_data/data_polynya_40/'
    # fixed_num = 40
    # polynya_detecting_fixed(datapd, p_polynya_save_fix, p_ice, area_threshold, fixed_num)

    print()
    print('Start sea ice thickness data preprocessing ...')
    save_new_dataset_sithickness(datapd, p_sith, p_nc, selected_month, southlat, newx=newx)
    
if __name__ == "__main__":
    main()