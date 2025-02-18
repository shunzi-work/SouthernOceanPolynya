###############
# 
# # Python script for CMIP6 hfds data preprocessing 
# 
###############

from myfunctions import *
from DataPre_siconc import get_new_dataset, save_new_dataset
import gc

def main():
    # filter some warning messages
    import warnings
    warnings.filterwarnings("ignore")
    
    datapd = pd.read_csv('List_model.csv')
    p_hfds = '../../SO_data/data_hfds/'
    
    
    p_nc = '../../../data/model/CMIP6/'
    selected_month = 9
    southlat = -40 

    print('Start mlotst data preprocessing ...')
    save_new_dataset(datapd, p_hfds, p_nc, selected_month, southlat, 'hfds')
    

if __name__ == "__main__":
    main()