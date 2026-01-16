###############
# 
# # Python script for CMIP6 so, thetao, mlotst data 
# # preprocessing and mld calculation
# 
###############

from myfunctions import *
from IceDataPre import get_new_dataset, save_new_dataset
import gc
from SeaDataPre import fix_time_coords

def exclude_bottom_error(da, levname, dep = 6700):
    return da.sel({levname: slice(0, dep)})


def cal_mld_from_ts_data(data_info, p_nc, selected_month, southlat):
    dst = get_new_dataset(data_info, p_nc, selected_month, southlat, 'thetao')
    if isinstance(dst, xr.Dataset):
        dss = get_new_dataset(data_info, p_nc, selected_month, southlat, 'so')
        if isinstance(dss, xr.Dataset):
            levname = data_info['zname']
            dat = check_lev_unit(levname, dst.thetao)
            das = check_lev_unit(levname, dss.so)

            dat = exclude_bottom_error(dat, levname)
            das = exclude_bottom_error(das, levname)
            
            da_sigma0 = gsw.sigma0(das, dat)
            try:
                da_mld = cal_mld(da_sigma0, levname)
                mld_new = xr.Dataset(
                    data_vars = {
                        'mld': da_mld,
                        'areacello':dst.areacello,  
                        'newlat': dst.newlat,     
                        'newlon': dst.newlon,
                    }
                )
                return mld_new
            except Exception as e:
                print(e)
                return None
        else:
            return None
    else:
        return None

def save_mld_calculation(datapd, p_save, p_nc, selected_month, southlat):
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        print("{} {}".format(i, name), end = '...')
        if ispickleexists(name, p_save):
            print("[o] data exist.")
            continue
        
        mld = cal_mld_from_ts_data(datapd.iloc[i], p_nc, selected_month, southlat)
        
        if isinstance(mld, xr.Dataset):
            mld = mld.load()
            savepickle(name, p_save, mld)
            print("[*] Saved.")   
            gc.collect()

def main():
    # filter some warning messages
    import warnings
    warnings.filterwarnings("ignore")
    
    datapd = pd.read_csv('List_model.csv')
    p_mlotst = '../../SO_data/data_mlotst/'
    p_mld = '../../SO_data/data_mld/'
    
    p_nc = '../../data/CMIP6/'
    selected_month = 9
    southlat = -40 

    print('Start mlotst data preprocessing ...')
    save_new_dataset(datapd, p_mlotst, p_nc, selected_month, southlat, 'mlotst')
    print('Finish mlotst data preprocessing.')
    print()
    print('Start mld calculation ...')
    save_mld_calculation(datapd, p_mld, p_nc, selected_month, southlat)
    
    p_ice = '../../SO_data/data_siconc_w_area/'
    for name in ['CESM2', 'CESM2-WACCM-FV2']:
        if ispickleexists(name, p_mld):
            fix_time_coords(name, p_mld, p_ice)
        if ispickleexists(name, p_mlotst):
            fix_time_coords(name, p_mlotst, p_ice)
    

if __name__ == "__main__":
    main()