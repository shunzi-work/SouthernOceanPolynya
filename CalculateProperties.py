###############
# 
# # Python script for calculating some properties 
# # reagarding polynya and convection each year
# 
###############

from myfunctions import *

def calculate_polynya_area(dspolynya, dsice):
    polynya_area_ann = dsice.areacello.where(dspolynya>0).sum((dsice.areacello.dims[0],dsice.areacello.dims[1]))/1e12
    freq = polynya_area_ann.where(polynya_area_ann>0).count()/len(polynya_area_ann.time)
    polynya_area_total = dsice.areacello.where(dspolynya.mean('time')>0).sum()/1e12
    ice_mean_not0 = (dsice.siconc.where(dsice.siconc.max("time")>0)*dsice.areacello).sum()/dsice.areacello.where(dsice.siconc.max("time")>0).sum()/len(dsice.time)
    return [ice_mean_not0.item(), polynya_area_ann.mean().item(), polynya_area_ann.max().item(), polynya_area_total.item(), freq.item()]


def check_timerange(dsice, ds):
    if len(dsice.time) != len(ds.time):
        return False
    elif (dsice.time != ds.time).any():
        return False
    else:
        return True

def open_mld(p_mlotst, p_mld, name):
    if ispickleexists(name, p_mlotst):
        ds = openpickle(name, p_mlotst)
        ds = change_start_x(ds, 135)
        da = ds.mlotst
    else:
        ds = openpickle(name, p_mld)
        ds = change_start_x(ds, 135)
        da = ds.mld
    if len(da.time)>500:
        da = da.isel(time = slice(-500, None))
    return da, ds

def calculate_convection_area(damld, daarea, dsice):
    if check_timerange(dsice, damld):
        convection_area_ann = daarea.where(damld>=2000).sum((daarea.dims[0], daarea.dims[1]))/1e12
        freq = convection_area_ann.where(convection_area_ann>0).count()/len(convection_area_ann)
        convection_area_total = daarea.where(damld.where(damld>=2000).mean('time')>0).sum()/1e12
        return [convection_area_ann.mean().item(), convection_area_ann.max().item(), convection_area_total.item(), freq.item()]
    else:
        print("time range different.")
        return [None, None, None, None]

def open_tos(name, p_tos):
    return open_multivarname(name, p_tos, ['thetao', 'tos'])

def open_sos(name, p_sos):
    return open_multivarname(name, p_sos, ['so', 'sos'])

def open_thick(name, p_thick):
    return open_multivarname(name, p_thick, ['sithick', 'sivol'])
    
def open_multivarname(mname, path_data, varnames):
    ds = openpickle(mname, path_data)
    ds = change_start_x(ds, 135)
    for varn in varnames:
        if varn in ds:
            da = ds[varn]
            break
    if len(da.time)>500:
        da = da.isel(time = slice(-500, None))
    return da, ds

def open_ann_stdata(mname, path_data, varnames):
    ds = openpickle(mname, path_data)
    ds = change_start_x(ds, 135)
    for varn in varnames:
        if varn in ds:
            da = ds[varn]
            break
    daann = da.groupby("time.year").mean("time")
    danew = daann.rename({'year': 'time'})
    danew['time'] = select_month(da, 9).time
    if len(danew.time)>500:
        danew = danew.isel(time = slice(-500, None))
    return danew

def cal_properties(da, newlat, dapolynya, damld):
    mean55 = da.where(newlat<=-55).mean()
    p_area = dapolynya.mean('time')
    c_area = damld.where(damld >= 2000).mean('time')
    if len(da[da.dims[-1]]) != len(dapolynya[dapolynya.dims[-1]]):
        print("ice coords not match.", end = '...')
        mean_p = None
    else:
        try:
            mean_po = da.where(p_area>0).mean()
            mean_p = mean_po.item()
        except Exception as E:
            print(E, end = '...')
            mean_p = None
    if len(da[da.dims[-1]]) != len(damld[damld.dims[-1]]):
        print("ocean coords not match.", end = '...')
        mean_c = None
    else:
        try:
            mean_co = da.where(c_area>0).mean()
            mean_c = mean_co.item()
        except Exception as E:
            print(E, end = '...')
            mean_c = None
    return [mean55.item(), mean_p, mean_c]

def add_new_properties(da, ds, dspolynya, damld, ps):
    vals = cal_properties(da, ds, dspolynya, damld)
    for val in vals:
        ps.append(val)
    return ps

