###############
# 
# # Python script for CMIP6 so, thetao, mlotst data 
# # preprocessing and mld calculation
# 
###############

from myfunctions import *
from DataPre_siconc import get_new_dataset, save_new_dataset
import gc

def cal_mld_from_ts_data(data_info, p_nc, selected_month, southlat):
    dst = get_new_dataset(data_info, p_nc, selected_month, southlat, 'thetao')
    if dst:        
        dss = get_new_dataset(data_info, p_nc, selected_month, southlat, 'so')
        if dss:
            levname = data_info['zname']
            dat = check_lev_unit(levname, dst.thetao)
            das = check_lev_unit(levname, dss.so)
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
    

    print('Finish siconc data preprocessing.')
    print()
    print('Start mld calculation ...')
    save_mld_calculation(datapd, p_mld, p_nc, selected_month, southlat)
    

if __name__ == "__main__":
    main()