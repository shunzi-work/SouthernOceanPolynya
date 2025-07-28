# ##############
#
# # Python script for calculating some properties 
# # reagarding polynya and convection each year
#
# ##############

from myfunctions import *
import gc

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


def main():
    p_polynya = '../../SO_data/data_polynya_mean/'
    p_ice = '../../SO_data/data_siconc_w_area/'
    p_mlotst = '../../SO_data/data_mlotst/'
    p_mld = '../../SO_data/data_mld/'
    p_tos = '../../SO_data/data_sst/'
    p_sos = '../../SO_data/data_sos/'
    p_thick = '../../SO_data/data_thick/'
    p_hfds = '../../SO_data/data_hfds/'
    p_tann = '../../SO_data/data_sst_ann/'
    p_sann = '../../SO_data/data_sos_ann/'
    
    mycols = [
        'name', 'p_criteria', \
        'p_mean', 'p_max', 'p_total', 'p_freq', \
        'c_mean', 'c_max', 'c_total', 'c_freq', \
        't_mean55', 't_mean_p', 't_mean_c', \
        's_mean55', 's_mean_p', 's_mean_c', \
        'ice_mean55', 'ice_mean_p', 'ice_mean_c', \
        'thick_m55', 'thick_mp', 'thick_mc', \
        'hfds_m55', 'hfds_mp', 'hfds_mc', \
        'mld_m55', 'mld_mp', 'mld_mc', \
        'tann_55', 'tann_p', 'tann_c', \
        'sann_55', 'sann_p', 'sann_c' 
    ]
    
    datapd = pd.read_csv('List_model.csv')
    df = pd.DataFrame(columns=mycols)
    df['name'] = datapd['source_id']
    for i in range(0, len(df)):
        ps = []
        name = df.at[i, 'name']
        print(name, end = '...')
        dsice = openpickle(name, p_ice)
        daice = dsice.siconc
        dspolynya = openpickle(name, p_polynya)
        polynya_ps = calculate_polynya_area(dspolynya, dsice)
        for val in polynya_ps:
            ps.append(val)
        damld, dsmld = open_mld(p_mlotst, p_mld, name)
        conv_ps = calculate_convection_area(damld, dsmld.areacello, dsice)
        for val in conv_ps:
            ps.append(val)
        dasst, dssst = open_tos(name, p_tos)
        dasss, dssss = open_sos(name, p_sos)
        dat_ann = open_ann_stdata(name, p_tann, ['tos', 'thetao'])
        das_ann = open_ann_stdata(name, p_sann, ['sos', 'so'])
        
        if dasst.dims[-1] != dspolynya.dims[-1]:
            print(dasst.dims, dspolynya.dims, end = '...')
            if len(dasst[dasst.dims[-1]]) == len(dspolynya[dspolynya.dims[-1]]):
                dspolynya = copy_xy(damld, dspolynya)
                daice = copy_xy(damld, daice)
        print('sst', end = '...')
        ps = add_new_properties(dasst, dssst.newlat, dspolynya, damld, ps)
        print('sss', end = '...')
        ps = add_new_properties(dasss, dssss.newlat, dspolynya, damld, ps)
        print('ice', end = '...')
        if len(dasst[dasst.dims[-1]]) != len(dspolynya[dspolynya.dims[-1]]):
            ps = add_new_properties(daice, dsice.newlat, dspolynya, damld, ps)
        else:
            ps = add_new_properties(daice, dsmld.newlat, dspolynya, damld, ps)
        print('thick', end = '...')
        if ispickleexists(name, p_thick):
            dathick, dsthick = open_thick(name, p_thick)
            dnewlat = dsmld.newlat
            if dathick.dims[-1] != damld.dims[-1]:
                print(dathick.dims, damld.dims, end = '...')
                if len(dathick[dathick.dims[-1]]) == len(damld[damld.dims[-1]]):
                    dathick = copy_xy(damld, dathick)
                else:
                    dnewlat = dsthick.newlat
            ps = add_new_properties(dathick, dnewlat, dspolynya, damld, ps)
        else:
            for k in range(3):
                ps.append(None)
        print('hfds', end = '...')
        if ispickleexists(name, p_hfds):
            dshfds = openpickle(name, p_hfds)
            dshfds = change_start_x(dshfds, 135)
            ps = add_new_properties(dshfds.hfds, dshfds.newlat, dspolynya, damld, ps)
        else:
            for k in range(3):
                ps.append(None)
        print('mld', end = '...')
        ps = add_new_properties(damld, dsmld.newlat, dspolynya, damld, ps)
        print('tos ann', end = '...')
        ps = add_new_properties(dat_ann, dssst.newlat, dspolynya, damld, ps)
        print('sos_ann', end = '...')
        ps = add_new_properties(das_ann, dssss.newlat, dspolynya, damld, ps)

        df.iloc[i, 1:] = ps
        print('')
    df.to_csv('properties.csv', index=False)

if __name__ == "__main__":
    main()
