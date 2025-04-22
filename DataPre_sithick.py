###############
# 
# # Python script for CMIP6 sea ice thickness data preprocessing
# 
###############

from myfunctions import *
from DataPre_siconc import get_new_dataset
import gc

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
    p_sith = '../../SO_data/data_thick/'
    newx = 135
    
    p_nc = '../../data/CMIP6/'
    selected_month = 9
    southlat = -40 

    print('Start sea ice data preprocessing ...')
    save_new_dataset_sithickness(datapd, p_sith, p_nc, selected_month, southlat, newx=newx)

    

if __name__ == "__main__":
    main()