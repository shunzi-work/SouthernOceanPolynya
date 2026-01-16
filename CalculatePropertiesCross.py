# ##############
#
# # Python script for calculating some properties
# # reagarding polynya and convection each year
#
# ##############

from myfunctions import *
from CalculateProperties import *
import gc

def calculate_cp_cross_rate(dspolynya, dsice, damld, dsmld):
    polynya_area_ann = dsice.areacello.where(dspolynya>0).sum((dsice.areacello.dims[0],dsice.areacello.dims[1]))/1e12
    convection_area_ann = dsmld.areacello.where(damld>=2000).sum((dsmld.areacello.dims[0], dsmld.areacello.dims[1]))/1e12
    if len(damld[damld.dims[-1]]) != len(dspolynya[dspolynya.dims[-1]]):
        print("[!] coords not match.")
        cross_area_ann = None
    else:
        try:
            cross_area_ann = dsmld.areacello.where(dspolynya>0).where(damld>=2000).sum((dsmld.areacello.dims[0], dsmld.areacello.dims[1]))/1e12
        except Exception as E:
            print(E, end = '...')
            cross_area_ann = None
        
    if isinstance(cross_area_ann, xr.DataArray):
        polynya_area_ann_p = polynya_area_ann.where(polynya_area_ann>0, drop = True)
        cross_area_ann_p = cross_area_ann.sel(time = polynya_area_ann_p.time)
        p_cross_rate = cross_area_ann_p/polynya_area_ann_p
        conv_area_ann_p = convection_area_ann.where(convection_area_ann>0, drop = True)
        cross_area_ann_c = cross_area_ann.sel(time = conv_area_ann_p.time)
        c_cross_rate = cross_area_ann_c/conv_area_ann_p
        return polynya_area_ann, convection_area_ann, p_cross_rate.mean().item(), c_cross_rate.mean().item()
    else:
        return polynya_area_ann, convection_area_ann, None, None


def main():
    import warnings
    warnings.filterwarnings("ignore")
    
    data_save_path= '../../SO_data/data_within_model/'
    p_polynya = '../../SO_data/data_polynya_mean/'
    p_ice = '../../SO_data/data_siconc_w_area/'
    p_mlotst = '../../SO_data/data_mlotst/'
    p_mld = '../../SO_data/data_mld/'

    datapd = pd.read_csv('List_model.csv')
    df = pd.DataFrame(columns=["name", "p_cross_rate", "c_cross_rate", "p_events", "c_events", "p_cross_c_events"])
    df["name"] = datapd['source_id']

    for i in range(0, len(datapd)):
        gc.collect()
        name = datapd.at[i, 'source_id']
        print(">>> {} {}".format(i, name), end = '...')
        
        dsice = openpickle(name, p_ice)
        daice = dsice.siconc
        dspolynya = openpickle(name, p_polynya)
        damld, dsmld = open_mld(p_mlotst, p_mld, name)

        if name == 'CESM2-WACCM-FV2':
            damld['time'] = dspolynya.time
            dsmld['time'] = dspolynya.time
        
        # if name == 'NorESM2-LM':
        #     print(damld.time)
        #     print(dspolynya.time)
        #     break

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

        polynya_area_ann, convection_area_ann, p_cross_rate, c_cross_rate = calculate_cp_cross_rate(dspolynya, dsice, damld, dsmld)
        df["p_cross_rate"].at[i] = p_cross_rate
        df["c_cross_rate"].at[i] = c_cross_rate

        p_events = polynya_area_ann.where(polynya_area_ann > 0.1)
        c_events = convection_area_ann.where(convection_area_ann > 0.1)

        p_cross_c_events = p_events.where(c_events > 0).count().item()/len(p_events.time)
        df["p_events"].at[i] = p_events.count().item()/len(p_events.time)
        df["c_events"].at[i] = c_events.count().item()/len(c_events.time)
        df["p_cross_c_events"].at[i] = p_cross_c_events

        print('')

    df.to_csv('CrossRate.csv', index=False)


if __name__ == "__main__":
    main()
