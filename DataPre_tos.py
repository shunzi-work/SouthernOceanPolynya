###############
# 
# # Python script for CMIP6 tos data preprocessing
# 
###############

from myfunctions import *
from DataPre_siconc import get_new_dataset, save_new_dataset
from DataPre_sos import save_new_dataset_surf
import gc

def main():
    # filter some warning messages
    import warnings
    warnings.filterwarnings("ignore")
    
    datapd = pd.read_csv('List_model.csv')
    p_sst = '../../SO_data/data_sst/'
    newx = 135
    
    p_nc = '../../../data/model/CMIP6/'
    selected_month = 9
    southlat = -40 

    print('Start tos data preprocessing ...')
    save_new_dataset(datapd, p_sst, p_nc, selected_month, southlat, 'tos')

    print('')
    print('Start thetao data preprocessing and save only surface ...')
    save_new_dataset_surf(datapd, p_sst, p_nc, selected_month, southlat, 'thetao', 'tos', newx=False)


if __name__ == "__main__":
    main() 


