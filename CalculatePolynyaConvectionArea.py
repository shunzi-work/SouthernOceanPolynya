###############
# 
# # Python script for calculating polynya and
# # convection area each year
# 
###############

from myfunctions import *
import gc

def calculate_polynya_area(p_polynya, p_ice, name):
    dsice = openpickle(name, p_ice)
    dspolynya = openpickle(name, p_polynya)
    polynya_area_ann = dsice.areacello.where(dspolynya>0).sum((dsice.areacello.dims[0],dsice.areacello.dims[1]))/1e12
    freq = polynya_area_ann.where(polynya_area_ann>0).count()/len(polynya_area_ann)
    polynya_area_total = dsice.areacello.where(dspolynya.mean('time')>0).sum()/1e12
    if name == 'KIOST-ESM':
        print(len(polynya_area_ann))
        

    return [polynya_area_ann.mean().item(), polynya_area_ann.max().item(), polynya_area_total.item(), freq.item()]

def calculate_convection_area(p_mlotst, p_mld, name):
    if ispickleexists(name, p_mlotst):
        ds = openpickle(name, p_mlotst)
    else:
        ds = openpickle(name, p_mld)
    if 'mlotst' in ds:
        damld = ds.mlotst
    else:
        damld = ds.mld
    convection_area_ann = ds.areacello.where(damld>=2000).sum((ds.areacello.dims[0], ds.areacello.dims[1]))/1e12
    freq = convection_area_ann.where(convection_area_ann>0).count()/len(convection_area_ann)
    convection_area_total = ds.areacello.where(damld.where(damld>=2000).mean('time')>0).sum()/1e12
    if name == 'KIOST-ESM':
        print(len(convection_area_ann))
    return [convection_area_ann.mean().item(), convection_area_ann.max().item(), convection_area_total.item(), freq.item()]
    
def main():
    p_polynya = '../../SO_data/data_polynya/'
    p_ice = '../../SO_data/data_siconc_w_area/'
    p_mlotst = '../../SO_data/data_mlotst/'
    p_mld = '../../SO_data/data_mld/'
    datapd = pd.read_csv('List_model.csv')
    area_list = []
    for i in range(0, len(datapd)):
        name = datapd.at[i, 'source_id']
        p_area = calculate_polynya_area(p_polynya, p_ice, name)
        c_area = calculate_convection_area(p_mlotst, p_mld, name)
        area_list.append([name] + p_area + c_area)
    area_df = pd.DataFrame(area_list, columns=['name', 'p_mean', 'p_max', 'p_total', 'p_freq', 'c_mean', 'c_max', 'c_total', 'c_freq'])
    # print(area_df)
    # area_df.to_csv('polynya_convection_areas.csv', index=False)

if __name__ == "__main__":
    main()
