###############
# 
# # Python script for CMIP6 sos and tos
# # data preprocessing
# 
###############

from myfunctions import *
from DataPre_siconc import get_new_dataset, save_new_dataset
import gc

def read_surface_data(datainfo, p_nc, selected_month, southlat, dataname, newx):
    new = get_new_dataset(datainfo, p_nc, selected_month, southlat, dataname, newx=newx)
    newds_surf = new.isel({datainfo['zname']:0})
    ds = newds_surf.reset_coords(datainfo['zname'], drop = True)
    return ds

def save_new_dataset_surf(datapd, p_save, p_nc, selected_month, southlat, dataname, newx=False):
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        print("{} {}".format(i, name), end = '...')
        if ispickleexists(name, p_save):
            print("[o] data exist.")
            continue
        new_ds_south = read_surface_data(datapd.iloc[i], p_nc, selected_month, southlat, dataname, newx=newx)
        if isinstance(new_ds_south, xr.Dataset):
            if (name in ['SAM0-UNICON', 'CAS-ESM2-0']) and (dataname == 'siconc'):
                # SAM0-UNICON, CAS-ESM2-0 with one year less in mld data than in ice data
                new_ds_south = new_ds_south.isel(time = slice(0, -1))
            new_ds_south = new_ds_south.load()
            savepickle(name, p_save, new_ds_south)
            print("[*] Saved.")
            gc.collect()
        tempfile_path = p_nc + '*' + '_temp.nc'
        matching_files_temp = glob.glob(tempfile_path)
        if len(matching_files_temp)>0:
            for mf in matching_files_temp:
                os.remove(mf)

def main():
    # filter some warning messages
    import warnings
    warnings.filterwarnings("ignore")
    
    datapd = pd.read_csv('List_model.csv')
    p_sos = '../../SO_data/data_sos/'
    p_sos_ann = '../../SO_data/data_sos_ann/'
    newx = 135
    
    p_nc = '../../data/CMIP6/'
    selected_month = 9
    southlat = -40 

    print('Start sos data preprocessing ...')
    save_new_dataset(datapd, p_sos, p_nc, selected_month, southlat, 'sos', newx)
    save_new_dataset(datapd, p_sos_ann, p_nc, 0, southlat, 'sos', newx)

    print('')
    print('Start so data preprocessing and save only surface ...')
    save_new_dataset_surf(datapd, p_sos, p_nc, selected_month, southlat, 'so', newx)
    save_new_dataset_surf(datapd, p_sos_ann, p_nc, 0, southlat, 'so', newx)
    
    print('')
    p_sst = '../../SO_data/data_sst/'
    p_sst_ann = '../../SO_data/data_sst_ann/'
    print('Start tos data preprocessing ...')
    save_new_dataset(datapd, p_sst, p_nc, selected_month, southlat, 'tos')
    save_new_dataset(datapd, p_sst_ann, p_nc, 0, southlat, 'tos')
    
    print('')
    print('Start thetao data preprocessing and save only surface ...')
    save_new_dataset_surf(datapd, p_sst, p_nc, selected_month, southlat, 'thetao', newx)
    save_new_dataset_surf(datapd, p_sst_ann, p_nc, 0, southlat, 'thetao', newx)
    
if __name__ == "__main__":
    main()

