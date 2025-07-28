# ##############
#
# # Python script for calculating some properties 
# # reagarding polynya and convection each year
#
# ##############

from myfunctions import *
from CalculateProperties import *
import gc

# def calculate_cp_cross_rate(dspolynya, dsice, damld, daarea):
#     polynya_area_ann = dsice.areacello.where(dspolynya>0).sum((dsice.areacello.dims[0],dsice.areacello.dims[1]))/1e12
#     convection_area_ann = daarea.where(damld>=2000).sum((daarea.dims[0], daarea.dims[1]))/1e12
#     if len(damld[damld.dims[-1]]) != len(dspolynya[dspolynya.dims[-1]]):
#         print("[!] ice coords not match.")
#         cross_area_ann = None
#     else:
#         try:
#             cross_area_ann = daarea.where(dspolynya>0).where(damld>=2000).sum((daarea.dims[0], daarea.dims[1]))/1e12
#         except Exception as E:
#             print(E, end = '...')
#             cross_area_ann = None
#     if isinstance(cross_area_ann, xr.DataArray):
#         polynya_area_ann_p = polynya_area_ann.where(polynya_area_ann>0, drop = True)
#         cross_area_ann_p = cross_area_ann.sel(time = polynya_area_ann_p.time)
#         p_cross_rate = cross_area_ann_p/polynya_area_ann_p
#         conv_area_ann_p = convection_area_ann.where(convection_area_ann>0, drop = True)
#         cross_area_ann_c = cross_area_ann.sel(time = conv_area_ann_p.time)
#         c_cross_rate = cross_area_ann_c/conv_area_ann_p
#         return [polynya_area_ann, convection_area_ann, p_cross_rate.mean().item(), c_cross_rate.mean().item()]
#     else:
#         return [polynya_area_ann, convection_area_ann, None, None]
    
# def calculate_var_corr(daregion, var1, var2):
#     var1_region = get_region(var1, daregion, 'region')
#     var2_region = get_region(var2, daregion, 'region')
#     if isinstance(var1_region, xr.DataArray) and isinstance(var2_region, xr.DataArray):
#         return np.corrcoef(var1_region, var2_region)[0, 1]
#     else:
#         return None

def get_region(da, region, region_name):
    if not isinstance(da, xr.DataArray):
        da_region = None
    elif len(da[da.dims[-1]]) != len(region[region.dims[-1]]):
        print("[!] {} coords not match.".format(region_name))
        da_region = None
    else:
        try:
            da_region = da.where(region>0).mean((region.dims[1], region.dims[0]))
        except Exception as E:
            print(E, end = '...')
            da_region = None
        return da_region

def get_region_diff(da, newlat, region, region_name):
    if not isinstance(da, xr.DataArray):
        da_region = None
    elif len(da[da.dims[-1]]) != len(region[region.dims[-1]]):
        print("[!] {} coords not match.".format(region_name))
        da_region = None
    else:
        try:
            mean55 = da.where(newlat<=-55).mean((region.dims[1], region.dims[0]))
            da_anom = da - mean55
            da_region = da_anom.where(region>0).mean((region.dims[1], region.dims[0])) 
        except Exception as E:
            print(E, end = '...')
            da_region = None
        return da_region

# def cal_properties_corr(da1, da2):
#     if isinstance(da1, xr.DataArray) and isinstance(da2, xr.DataArray):
#         if check_timerange(da1, da2):
#             da1std = np.std(da1)
#             da2std = np.std(da2)
#             if (da1std == 0) or (da2std == 0):
#                 cc = None
#             else:
#                 cc = np.corrcoef(da1, da2)[0, 1]
#         else:
#             print("[!] time range different no cc.")
#             cc = None
#     else:
#         print("[!] One of the data is not a DataArray.")
#         cc = None
#     return cc

# def creat_vardict(var_names):
#     var_dict = {}
#     for v1 in range(0, len(var_names)):
#         for v2 in range(0, len(var_names)):
#             if v1 < v2:
#                 newcolname1 = var_names[v1] + '_vs_' + var_names[v2] + '_c'
#                 newcolname2 = var_names[v1] + '_vs_' + var_names[v2] + '_p'
#                 var_dict[newcolname1] = [var_names[v1], var_names[v2], 'c']
#                 var_dict[newcolname2] = [var_names[v1], var_names[v2], 'p']
#     return var_dict


def readvardata(var_name, name, data_path_pre, damld):
    # var_names = ['carea', 'parea', 'sic', 'mld', 'thick', 'hfds', 'sst', 'sss', 'sann', 'tann']
    data_path = data_path_pre + var_name + '/'
    dataarray_names = ['thetao', 'tos', 'sos', 'so', 'hfds', 'sithick', 'sivol']
    if ispickleexists(name, data_path):
        if var_name == 'sst_ann' or var_name == 'sos_ann':
            da = open_ann_stdata(name, data_path, dataarray_names)
        else:
            da, __ = open_multivarname(name, data_path, dataarray_names)
    else:
        print("[!] No data in {}".format(data_path))
        da = None
        
    if isinstance(da, xr.DataArray):
        if var_name == 'thick':
            if da.dims[-1] != damld.dims[-1]:
                print(da.dims, damld.dims, end = '...')
                if len(da[da.dims[-1]]) == len(damld[damld.dims[-1]]):
                    da = copy_xy(damld, da)
        if (var_name == 'hfds') and (name == 'SAM0-UNICON'):
            da = da.isel(time = slice(-500, -1))
    return da

# def get_vardata_region(varname, name, datadict, datadict2, area_data, area_type, data_path_pre, damld):
#     if varname in ['carea', 'parea']:
#         var_region = datadict[varname]
#     else:
#         if varname in datadict2:
#             var_region = datadict2[varname]
#         else:
#             if varname in ['sic', 'mld']:
#                 var = datadict[varname]
#             else:
#                 var = readvardata(varname, name, data_path_pre, damld)
#             var_region = get_region(var, area_data, area_type)
#     return var_region

# def update_list_and_dict(var_region, varname, empty_var_list, datadict):
#     if not isinstance(var_region, xr.DataArray):
#         empty_var_list.append(varname)
#     else:
#         datadict[varname] = var_region
#     return empty_var_list, datadict

# def get_vardata_region_from_colname(varnames, varname_dict, name, datadict, datadict2, area_data, area_type, data_path_pre, damld, empty_var_list):
#     var1_region = get_vardata_region(varname_dict[varnames][0], name, datadict, datadict2, area_data, area_type, data_path_pre, damld)
#     var2_region = get_vardata_region(varname_dict[varnames][1], name, datadict, datadict2, area_data, area_type, data_path_pre, damld)

#     empty_var_list, datadict2 = update_list_and_dict(var1_region, varname_dict[varnames][0], empty_var_list, datadict2)
#     empty_var_list, datadict2 = update_list_and_dict(var2_region, varname_dict[varnames][1], empty_var_list, datadict2)
#     return var1_region, var2_region, empty_var_list, datadict2
    

def main():
    import warnings
    warnings.filterwarnings("ignore")
    
    data_save_path= '../../SO_data/data_within_model/'
    data_path_pre = '../../SO_data/data_'
    p_polynya = '../../SO_data/data_polynya_mean/'
    p_ice = '../../SO_data/data_siconc_w_area/'
    p_mlotst = '../../SO_data/data_mlotst/'
    p_mld = '../../SO_data/data_mld/'

    var_names = ['thick', 'hfds', 'sst', 'sos', 'sst_ann', 'sos_ann']
    datapd = pd.read_csv('List_model.csv')
    for i in range(0, len(datapd)):
        gc.collect()
        name = datapd.at[i, 'source_id']
        print(">>> {} {}".format(i, name), end = '...')

        if ispickleexists(name, data_save_path):
            print("done.")
            continue
        
        dsice = openpickle(name, p_ice)
        daice = dsice.siconc
        dspolynya = openpickle(name, p_polynya)
        damld, dsmld = open_mld(p_mlotst, p_mld, name)

        p_area = dspolynya.mean('time')
        c_area = damld.where(damld >= 2000).mean('time')

        if p_area.sum().item() == 0.0:
            p_events = False
        else:
            p_events = True
        if c_area.sum().item() == 0.0:
            c_events = False
        else:
            c_events = True

        if not p_events and not c_events:
            print("no c and no p.")
            continue

        
        if (damld.dims[-1] != dspolynya.dims[-1]) or (damld.dims[-2] != dspolynya.dims[-2]):
            print(damld.dims, dspolynya.dims, end = '...')
            if len(damld[damld.dims[-1]]) == len(dspolynya[dspolynya.dims[-1]]):
                dspolynya = copy_xy(damld, dspolynya)
                daice = copy_xy(damld, daice)
                print("copy dims to ice grid")
            else:
                print("cords have different size")

        if not check_timerange(dsice, damld):
            print("[!] c and p time range different.")
            continue

        polynya_area_ann = dsice.areacello.where(dspolynya>0).sum((dsice.areacello.dims[0],dsice.areacello.dims[1]))/1e12
        convection_area_ann = dsmld.areacello.where(damld>=2000).sum((dsmld.areacello.dims[0], dsmld.areacello.dims[1]))/1e12
        
        ds = xr.Dataset({'polynya_area_ann': polynya_area_ann, 'convection_area_ann': convection_area_ann})
        print("sic", end = '...')
        ds['sic_convection'] = get_region(daice, c_area, "convection")
        ds['sic_polynya'] = get_region(daice, p_area, "polynya")
        
        ds['sic_diff_convection'] = get_region_diff(daice, dsice.newlat, c_area, "convection")
        ds['sic_diff_polynya'] = get_region_diff(daice, dsice.newlat, p_area, "polynya")

        print("mld", end = '...')
        ds['mld_convection'] = get_region(damld, c_area, "convection")
        ds['mld_diff_convection'] = get_region_diff(damld, dsmld.newlat, c_area, "convection")
        ds['mld_polynya'] = get_region(damld, p_area, "polynya")
        ds['mld_diff_polynya'] = get_region_diff(damld, dsmld.newlat, p_area, "polynya")

        for varname in var_names:
            print(varname, end='...')
            var = readvardata(varname, name, data_path_pre, damld)
            if isinstance(var, xr.DataArray):
                ds[varname+ '_polynya'] = get_region(var, p_area, "polynya")
                ds[varname + '_convection'] = get_region(var, c_area, "convection")
                if varname == 'thick':
                    ds[varname + '_diff_polynya'] = get_region_diff(var, dsice.newlat, p_area, "polynya")
                    ds[varname + '_diff_convection'] = get_region_diff(var, dsice.newlat, c_area, "convection")
                else:
                    ds[varname + '_diff_polynya'] = get_region_diff(var, dsmld.newlat, p_area, "polynya")
                    ds[varname + '_diff_convection'] = get_region_diff(var, dsmld.newlat, c_area, "convection")
            else:
                print("[!] {} data not found.".format(varname))
            del var
            gc.collect()

        for t in ds:
            if not isinstance(ds[t], xr.DataArray):
                ds.drop_vars(t)
        
        savepickle(name, data_save_path, ds)
        print('')
    

if __name__ == "__main__":
    main()
