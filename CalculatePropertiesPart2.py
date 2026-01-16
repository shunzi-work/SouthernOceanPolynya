# ##############
#
# # Python script for calculating some properties 
# # reagarding polynya and convection each year
#
# ##############

from myfunctions import *
from CalculateProperties import *
import gc

from scipy import stats


def cal_properties_diff(da, newlat, dapolynya, damld):
    mean55 = da.where(newlat<=-55).mean()
    p_area = dapolynya.mean('time')
    c_area = damld.where(damld >= 2000).mean('time')
    if len(da[da.dims[-1]]) != len(dapolynya[dapolynya.dims[-1]]):
        print("ice coords not match.", end = '...')
        mean_p_diff = None
    else:
        try:
            mean_p_diff = (da.where(p_area>0) - mean55).mean().item()
        except Exception as E:
            print(E, end = '...')
            mean_p_diff = None
    if len(da[da.dims[-1]]) != len(damld[damld.dims[-1]]):
        print("ocean coords not match.", end = '...')
        mean_c_diff = None
    else:
        try:
            mean_c_diff = (da.where(c_area>0) - mean55).mean().item()
        except Exception as E:
            print(E, end = '...')
            mean_c_diff = None
    return [mean_p_diff, mean_c_diff]

def add_new_properties_diff(da, ds, dspolynya, damld, ps):
    vals = cal_properties_diff(da, ds, dspolynya, damld)
    for val in vals:
        ps.append(val)
    return ps


def add_new_properties_pcdiff(da, dspolynya, damld, ps):
    if check_timerange(dspolynya, da):
        vals = calculate_pcdiff(dspolynya, damld, da)
    else:
        print('Time different')
        vals = [None, None]
    for val in vals:
        ps.append(val)
    return ps

def add_new_properties_pcdiff_ttest(da, dspolynya, damld, ps):
    if check_timerange(dspolynya, da):
        vals = calculate_pcdiff_ttest(dspolynya, damld, da)
    else:
        print('Time different')
        vals = [None, None]
    for val in vals:
        ps.append(val)
    return ps


def calculate_pcdiff(dspolynya, damld, da):
    if len(damld[damld.dims[-1]]) != len(dspolynya[dspolynya.dims[-1]]):
        print('different coords', end = '...')
        diff1 = None
        diff2 = None
        return [diff1, diff2]
    else:
        try:
            pandc = xr.where(dspolynya.where(damld >= 2000)>0, 1, 0)
            if np.isnan(pandc).all():
                diff1 = None
                diff2 = None
            else:
                nconvp = xr.where((dspolynya > 0) & (pandc == 0), 1, 0)
                npconv = xr.where((damld >= 2000) & (pandc == 0), 1, 0)
                
                pandc_da = da.where(pandc > 0).mean()
                nconvp_da = da.where(nconvp > 0).mean()
                npconv_da = da.where(npconv > 0).mean()
                
                diff1 = (pandc_da - nconvp_da).item()
                diff2 = (pandc_da - npconv_da).item()
            return [diff1, diff2]
        except Exception as E:
            print(E, end = '...')
            diff1 = None
            diff2 = None
        return [diff1, diff2]

def flatten_data(da):
    flattened = da.values.ravel()
    non_nan_data = flattened[~np.isnan(flattened)]
    return non_nan_data


def calculate_pcdiff_ttest(dspolynya, damld, da):
    if len(damld[damld.dims[-1]]) != len(dspolynya[dspolynya.dims[-1]]):
        print('different coords', end = '...')
        diff1 = None
        ttest = None
        return [diff1, ttest]
    else:
        try:
            pandc = xr.where(dspolynya.where(damld >= 2000)>0, 1, 0)
            if np.isnan(pandc).all():
                diff1 = None
                ttest = None
            else:
                nconvp = xr.where((dspolynya > 0) & (pandc == 0), 1, 0)
                
                pandc_da = da.where(pandc > 0)
                nconvp_da = da.where(nconvp > 0)

                diff1 = (pandc_da.mean() - nconvp_da.mean()).item()
                pandc_f = flatten_data(pandc_da)
                nconvp_f = flatten_data(nconvp_da)
                _, ttest = stats.ttest_ind(pandc_f, nconvp_f, equal_var=False)

            return [diff1, ttest]
        except Exception as E:
            print(E, end = '...')
            diff1 = None
            ttest = None
        return [diff1, ttest]



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
        'name', 't_mean_p_diff', 't_mean_c_diff', \
        's_mean_p_diff', 's_mean_c_diff', \
        'ice_mean_p_diff', 'ice_mean_c_diff', \
        'thick_mp_diff', 'thick_mc_diff', \
        'hfds_mp_diff', 'hfds_mc_diff', \
        'mld_mp_diff', 'mld_mc_diff', \
        'tann_p_diff', 'tann_c_diff', \
        'sann_p_diff', 'sann_c_diff' 
    ]

    mycols2 = [
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
        'sann_55', 'sann_p', 'sann_c', \
    ]

    mycols3 = [
        'name', 
        't_pcnc', 't_cpnp',\
        's_pcnc', 's_cpnp',\
        'ice_pcnc', 'ice_cpnp',\
        'thick_pcnc', 'thick_cpnp',\
        'hfds_pcnc', 'hfds_cpnp',\
        'mld_pcnc', 'mld_cpnp',\
        'tann_pcnc', 'tann_cpnp',\
        'sann_pcnc', 'sann_cpnp',\
    ]

    mycols4 = [
        'name', 
        't_pcnc', 't_test',\
        's_pcnc', 's_test',\
        'ice_pcnc', 'ice_test',\
        'thick_pcnc', 'thick_test',\
        'hfds_pcnc', 'hfds_test',\
        'mld_pcnc', 'mld_test',\
        'tann_pcnc', 'tann_test',\
        'sann_pcnc', 'sann_test',\
    ]
    
    datapd = pd.read_csv('List_model.csv')
    df = pd.DataFrame(columns=mycols)
    df['name'] = datapd['source_id']

    df2 = pd.DataFrame(columns=mycols2)
    df2['name'] = datapd['source_id']

    df3 = pd.DataFrame(columns=mycols3)
    df3['name'] = datapd['source_id']

    df4 = pd.DataFrame(columns=mycols4)
    df4['name'] = datapd['source_id']

    for i in range(0, len(df)):
        gc.collect()
        ps = []
        ps2 = []
        ps3 = []
        ps4 = []
        
        name = df.at[i, 'name']
        print(name, end = '...')
        dsice = openpickle(name, p_ice)
        daice = dsice.siconc
        dspolynya = openpickle(name, p_polynya)
            
        polynya_ps = calculate_polynya_area(dspolynya, dsice)
        for val in polynya_ps:
            ps2.append(val)

        damld, dsmld = open_mld(p_mlotst, p_mld, name)
        conv_ps = calculate_convection_area(damld, dsmld.areacello, dsice)
        for val in conv_ps:
            ps2.append(val)

        if not check_timerange(dsice, dsmld):
            if name == 'CESM2-WACCM-FV2':
                damld['time'] = dspolynya.time
                dsmld['time'] = dspolynya.time
            else:
                print('Time different')
                continue

        dasst, dssst = open_tos(name, p_tos)
        dasss, dssss = open_sos(name, p_sos)
        dat_ann = open_ann_stdata(name, p_tann, ['tos', 'thetao'])
        dat_ann['time'] = dasst['time']
        das_ann = open_ann_stdata(name, p_sann, ['sos', 'so'])
        das_ann['time'] = dasss['time']
        
        if dasst.dims[-1] != dspolynya.dims[-1]:
            print(dasst.dims, dspolynya.dims, end = '...')
            if len(dasst[dasst.dims[-1]]) == len(dspolynya[dspolynya.dims[-1]]):
                dspolynya = copy_xy(damld, dspolynya)
                daice = copy_xy(damld, daice)


        print('sst', end = '...')
        ps2 = add_new_properties(dasst, dssst.newlat, dspolynya, damld, ps2)
        ps = add_new_properties_diff(dasst, dssst.newlat, dspolynya, damld, ps)
        ps3 = add_new_properties_pcdiff(dasst, dspolynya, damld, ps3)
        ps4 = add_new_properties_pcdiff_ttest(dasst, dspolynya, damld, ps4)

        print('sss', end = '...')
        ps2 = add_new_properties(dasss, dssss.newlat, dspolynya, damld, ps2)
        ps = add_new_properties_diff(dasss, dssss.newlat, dspolynya, damld, ps)
        ps3 = add_new_properties_pcdiff(dasss, dspolynya, damld, ps3)
        ps4 = add_new_properties_pcdiff_ttest(dasss, dspolynya, damld, ps4)

        print('ice', end = '...')
        if len(dasst[dasst.dims[-1]]) != len(dspolynya[dspolynya.dims[-1]]):
            ps = add_new_properties_diff(daice, dsice.newlat, dspolynya, damld, ps)
            ps2 = add_new_properties(daice, dsice.newlat, dspolynya, damld, ps2)
            ps3 = add_new_properties_pcdiff(daice, dspolynya, damld, ps3)
            ps4 = add_new_properties_pcdiff_ttest(daice, dspolynya, damld, ps4)
        else:
            ps = add_new_properties_diff(daice, dsmld.newlat, dspolynya, damld, ps)
            ps2 = add_new_properties(daice, dsmld.newlat, dspolynya, damld, ps2)
            ps3 = add_new_properties_pcdiff(daice, dspolynya, damld, ps3)
            ps4 = add_new_properties_pcdiff_ttest(daice, dspolynya, damld, ps4)

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
            ps = add_new_properties_diff(dathick, dnewlat, dspolynya, damld, ps)
            ps2 = add_new_properties(dathick, dnewlat, dspolynya, damld, ps2)
            ps3 = add_new_properties_pcdiff(dathick, dspolynya, damld, ps3)
            ps4 = add_new_properties_pcdiff_ttest(dathick, dspolynya, damld, ps4)
        else:
            for k in range(3):
                if k < 2:
                    ps.append(None)
                    ps3.append(None)
                    ps4.append(None)
                ps2.append(None)

        print('hfds', end = '...')
        if ispickleexists(name, p_hfds):
            dshfds = openpickle(name, p_hfds)
            dshfds = change_start_x(dshfds, 135)
            if len(dshfds.time)>500:
                dshfds = dshfds.isel(time = slice(-500, None))
            ps = add_new_properties_diff(dshfds.hfds, dshfds.newlat, dspolynya, damld, ps)
            ps2 = add_new_properties(dshfds.hfds, dshfds.newlat, dspolynya, damld, ps2)
            ps3 = add_new_properties_pcdiff(dshfds.hfds, dspolynya, damld, ps3)
            ps4 = add_new_properties_pcdiff_ttest(dshfds.hfds, dspolynya, damld, ps4)
        else:
            for k in range(3):
                if k < 2:
                    ps.append(None)
                    ps3.append(None)
                    ps4.append(None)
                ps2.append(None)

        print('mld', end = '...')
        ps = add_new_properties_diff(damld, dsmld.newlat, dspolynya, damld, ps)
        ps2 = add_new_properties(damld, dsmld.newlat, dspolynya, damld, ps2)
        ps3 = add_new_properties_pcdiff(damld, dspolynya, damld, ps3)
        ps4 = add_new_properties_pcdiff_ttest(damld, dspolynya, damld, ps4)

        print('tos ann', end = '...')
        ps = add_new_properties_diff(dat_ann, dssst.newlat, dspolynya, damld, ps)
        ps2 = add_new_properties(dat_ann, dssst.newlat, dspolynya, damld, ps2)
        ps3 = add_new_properties_pcdiff(dat_ann, dspolynya, damld, ps3)
        ps4 = add_new_properties_pcdiff_ttest(dat_ann, dspolynya, damld, ps4)

        print('sos_ann', end = '...')
        ps = add_new_properties_diff(das_ann, dssss.newlat, dspolynya, damld, ps)
        ps2 = add_new_properties(das_ann, dssss.newlat, dspolynya, damld, ps2)
        ps3 = add_new_properties_pcdiff(das_ann, dspolynya, damld, ps3)
        ps4 = add_new_properties_pcdiff_ttest(das_ann, dspolynya, damld, ps4)

        df.iloc[i, 1:] = ps
        df2.iloc[i, 1:] = ps2
        df3.iloc[i, 1:] = ps3
        df4.iloc[i, 1:] = ps4
        print('')
    df.to_csv('properties_diff.csv', index=False)
    df2.to_csv('properties.csv', index=False)
    df3.to_csv('properties_pcdiff.csv', index=False)
    df4.to_csv('properties_pcdiff_ttest.csv', index=False)


def main2():
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
    
    mycols2 = [
        'name', 't_mean55', 't_std55', \
        'ice_mean55', 'ice_std55'
    ]
    
    datapd = pd.read_csv('List_model.csv')
    
    df2 = pd.DataFrame(columns=mycols2)
    df2['name'] = datapd['source_id']
    
    for i in range(0, len(df2)):
        gc.collect()
        ps2 = []
        name = df2.at[i, 'name']
        print(name, end = '...')
        dsice = openpickle(name, p_ice)
        daice = dsice.siconc
        dspolynya = openpickle(name, p_polynya)
    
        damld, dsmld = open_mld(p_mlotst, p_mld, name)
        conv_ps = calculate_convection_area(damld, dsmld.areacello, dsice)
    
        dasst, dssst = open_tos(name, p_tos)
        dasss, dssss = open_sos(name, p_sos)
        dat_ann = open_ann_stdata(name, p_tann, ['tos', 'thetao'])
        das_ann = open_ann_stdata(name, p_sann, ['sos', 'so'])
        
        # if dasst.dims[-1] != dspolynya.dims[-1]:
        #     print(dasst.dims, dspolynya.dims, end = '...')
        #     if len(dasst[dasst.dims[-1]]) == len(dspolynya[dspolynya.dims[-1]]):
        #         dspolynya = copy_xy(damld, dspolynya)
        #         daice = copy_xy(damld, daice)
    
        df2.loc[i, 't_mean55'] = dasst.where(dssst.newlat<=-55).mean().item()
        df2.loc[i, 'ice_mean55'] = daice.where(dsice.newlat<=-55).mean().item()
    
        df2.loc[i, 't_std55'] = dasst.where(dssst.newlat<=-55).std().item()
        df2.loc[i, 'ice_std55'] = daice.where(dsice.newlat<=-55).std().item()
    
    #     df.iloc[i, 1:] = ps
    #     df2.iloc[i, 1:] = ps2
        print('')
    # df.to_csv('properties_diff.csv', index=False)
    df2.to_csv('properties_sst_ice.csv', index=False)


if __name__ == "__main__":
    main()

