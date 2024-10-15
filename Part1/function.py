# Preprocessing for CMIP6 models
import warnings

import numpy as np
import xarray as xr
import pandas as pd

import os
import glob
import pickle

from pyproj import Geod

# def cmip6_renaming_dict():
#     """a universal renaming dict. Keys correspond to source id (model name)
#     and valuse are a dict of target name (key) and a list of variables that
#     should be renamed into the target."""

#     xyname1 = {
#         "xname" : "x",
#         "yname" : "y"
#         }
#     xyname2 = {
#         "xname" : "i",
#         "yname" : "j"
#         }
#     xyname3 = {
#         "xname" : "ni",
#         "yname" : "nj"
#         }
#     xyname4 = {
#         "xname" : "lon",
#         "yname" : "lat"
#         }
#     xyname5 = {
#         "xname" : "nlon",
#         "yname" : "nlat"
#     }

#     latlonname1 = {
#         "lonname" : "lon",
#         "latname" : "lat"
#     }
#     latlonname2 = {
#         "lonname" : "longitude",
#         "latname" : "latitude"
#     }
#     latlonname3 = {
#         "lonname" : "nav_lon",
#         "latname" : "nav_lat"
#     }

#     pair1 = {
#         "xyname" : xyname1, #x,y
#         "latlonname" : latlonname1, #lon, lat
#         "model" : {
#             "gn" : [
#                 "GFDL-CM4", 
#                 "GFDL-ESM4",
#                 "CNRM-CM6-1",
#                 "CNRM-CM6-1-HR",
#                 "CNRM-ESM2-1"
#             ],
#             "gr" : []
#         }
#     }

#     pair2 = {
#         "xyname" : xyname1, #x,y
#         "latlonname" : latlonname2, #longitude, latitude
#         "model" : {
#             "gn" : [
#                 "MIROC6",
#                 "MIROC-ES2L",
#                 "MRI-ESM2-0"
#             ],
#             "gr" : []
#         }
#     }

#     pair3 = {
#         "xyname" : xyname1, #x,y
#         "latlonname" : latlonname3, #nav_lon, nav_lat
#         "model" : {
#             "gn" : [
#                 "IPSL-CM5A2-INCA",
#                 "IPSL-CM6A-LR"
#             ],
#             "gr" : []
#         }
#     }

#     pair4 = {
#         "xyname" : xyname2,
#         "latlonname" : latlonname2,
#         "model" : {
#             "gn" : [
#                 "CAMS-CSM1-0",
#                 "ACCESS-ESM1-5",
#                 "ACCESS-CM2",
#                 "NESM3",
#                 "CanESM5",
#                 "CanESM5-1",
#                 "CanESM5-CanOE",
#                 "CMCC-CM2-SR5",
#                 "CMCC-ESM2",
#                 "EC-Earth3",
#                 "EC-Earth3-CC",
#                 "EC-Earth3-LR",
#                 "EC-Earth3-Veg",
#                 "EC-Earth3-Veg-LR",
#                 "HadGEM3-GC31-LL",
#                 "HadGEM3-GC31-MM",
#                 "UKESM1-0-LL",
#                 "FGOALS-g3",
#                 "CIESM",
#                 "TaiESM1",
#                 "SAM0-UNICON",
#                 "NorESM1-F",
#                 "NorESM2-MM",
#                 "NorCPM1",
#                 "MPI-ESM-1-2-HAM",
#                 "MPI-ESM1-2-HR",
#                 "MPI-ESM1-2-LR",
#                 ],
#             "gr" : []
#         }
#     }
    
#     pair5 = {
#         "xyname" : xyname4,
#         "latlonname" : latlonname2,
#         "model" : {
#             "gn" : [
#                 "BCC-CSM2-MR", 
#                 "BCC-ESM1"
#             ],
#             "gr" : []
#         }
#     }

#     pair6 = {
#         "xyname" : xyname4,
#         "latlonname" : None,
#         "model" : {
#             "gn" : [],
#             "gr" : [
#                 "GISS-E2-1-H",
#                 "GISS-E2-2-H",
#                 "E3SM-1-0",
#                 "E3SM-1-1",
#                 "E3SM-1-1-ECA",
#                 "E3SM-2-0"
#             ],
#             "gr1": [
#                 "INM-CM4-8",
#                 "INM-CM5-0"
#             ]
#         }
#     }

#     pair7 = {
#         "xyname" : xyname5,
#         "latlonname" : latlonname1,
#         "model" : {
#             "gn" : [
#                 "CESM2",
#                 "CESM2-FV2",
#                 "CESM2-WACCM",
#                 "CESM2-WACCM-FV2"
#             ],
#             "gr" : []
#         }
#     }
#     return rename_dict

def ispickleexists(n, p0):
    p = p0 + n + '.pickle'
    if os.path.exists(p):
        # print('    [o] {} exists.'.format(p))
        return True
    else:
        return False

def openpickle(n, p0):
    p = p0 + n + '.pickle'
    with open(p, 'rb') as df:
        d = pickle.load(df)
    return d

def savepickle(n, p0, sf):
    p = p0 + n + '.pickle'
    with open(p, 'wb') as wf:
        pickle.dump(sf, wf, pickle.HIGHEST_PROTOCOL)

def get_south(da, datapd, i, southlat = -40, newlatname =0):
    if pd.isna(datapd.at[i, 'latname']):
        da_south = da.sel({datapd.at[i, 'yname']:slice(-90, southlat)})
    else:
        latname = datapd.at[i, 'latname']
        if newlatname:
            latname = newlatname
        da_lat0 = da[latname].load()
        da_lat = da_lat0.where((da_lat0<=90) & (da_lat0>=-90))
        da_south = da.where(da_lat <= southlat, drop=True)
    return da_south

def get_sepsouth(mf, datapd, i, southlat = -40):
    ds = xr.open_mfdataset(mf, use_cftime=True)
    da = ds.siconc
    if 'type' in da.coords:
        da = da.reset_coords('type', drop = True)
    da_sep = da.isel(time=(da.time.dt.month == 9))
    da_south = get_south(da_sep, datapd, i, southlat)
    return da_south

def opennc(p0, datapd, i, southlat = -40):
    data_path = p0 + datapd.at[i, 'source_id'] + '_piControl_' + datapd.at[i, 'member_id'] + '*' + datapd.at[i, 'grid_label'] + '*' + '.nc'
    if datapd.at[i, 'source_id']  == 'NorESM2-LM' or datapd.at[i, 'source_id'] == 'NorESM2-MM':
        data_path = p0 + datapd.at[i, 'source_id'] + '_piControl_' + datapd.at[i, 'member_id'] + '*' + '.nc'
    matching_files = glob.glob(data_path)
    if len(matching_files)>200:    
        for t in range(len(matching_files)):
            da_south = get_sepsouth(matching_files[t], datapd, i, southlat = southlat)
            if t == 0:
                da_save = da_south.load()
            else:
                da_save0 = da_south.load()
                da_save = xr.concat([da_save, da_save0], dim="time")
        da_s = da_save
    else:
        da_south = get_sepsouth(matching_files, datapd, i, southlat = southlat)
        da_s = da_south.load()
    return da_s

def get_new_xy(d, datapd, j, newlonlat = 0):
    if newlonlat:
        newlon = d[newlonlat[0]]
        newlat = d[newlonlat[1]]
    else:
        if pd.isna(datapd.at[j, 'latname']):
            pltx0 = d[datapd.at[j, 'xname']]
            plty0 = d[datapd.at[j, 'yname']]
            newlon0, newlat0 = np.meshgrid(pltx0, plty0)
            newlon = xr.DataArray(newlon0, dims={datapd.at[j, 'yname']:plty0.values, datapd.at[j, 'xname']:pltx0.values})
            newlat = xr.DataArray(newlat0, dims={datapd.at[j, 'yname']:plty0.values, datapd.at[j, 'xname']:pltx0.values})
        else:
            newlon = d[datapd.at[j, 'lonname']]
            newlat = d[datapd.at[j, 'latname']]
    if len(np.shape(newlat)) > 2:
        newlon = newlon.isel(time = 0)
        newlat = newlat.isel(time = 0)
    return newlon, newlat

def copy_coords(icedata, areadata, datapd, i):
    newd = areadata.assign_coords(
        {
            datapd.at[i, 'latname']: icedata[datapd.at[i, 'latname']],
            datapd.at[i, 'lonname']: icedata[datapd.at[i, 'lonname']],
        }
    )
    return newd

def rename_xy(data_copyfrom, data_copyto):
    new = data_copyto.rename(
        {
            data_copyto.dims[len(data_copyto.dims)-1]:data_copyfrom.dims[len(data_copyfrom.dims)-1], 
            data_copyto.dims[len(data_copyto.dims)-2]:data_copyfrom.dims[len(data_copyfrom.dims)-2], 
        }
    )
    return new

def copy_coords_xy(data_copyfrom, data_copyto, changename = False):
    if changename:
        data_copyto = rename_xy(data_copyfrom, data_copyto)
    data_copyto[data_copyto.dims[len(data_copyto.dims)-1]] = data_copyfrom[data_copyfrom.dims[len(data_copyfrom.dims)-1]].values
    data_copyto[data_copyto.dims[len(data_copyto.dims)-2]] = data_copyfrom[data_copyfrom.dims[len(data_copyfrom.dims)-2]].values
    return data_copyto

def newxy_fmissingxy(dx, dy):
    dx = dx.where((dx>-361) & (dx<361))
    dy = dy.where((dy>-91) & (dy<91))
    newx0 = dx[~np.isnan(dx).any(axis=1)][0]
    newy0 = dy[:, ~np.isnan(dy).any(axis=0)][:,0]
    newx, newy = np.meshgrid(newx0, newy0)
    x = np.where(np.isnan(dx), newx, dx)
    y = np.where(np.isnan(dy), newy, dy)
    return x, y

def replace_missingxg(d, datapd, j, newlonlat = 0):
    if newlonlat:
        lon = d[newlonlat[0]]
        lat = d[newlonlat[1]]
    else:
        lon = d[datapd.at[j, 'lonname']]
        lat = d[datapd.at[j, 'latname']]
    if 'time' in lat.dims:
        lon = lon.isel(time = 0)
        lat = lat.isel(time = 0)
    newlon, newlat = newxy_fmissingxy(lon, lat)
    newd = d.assign_coords(
        {
            datapd.at[j, 'lonname']: (lon.dims, newlon),
            datapd.at[j, 'latname']: (lat.dims, newlat),
        }
    )
    return newd

def calculate_area_xy(icedata):
    g = Geod(ellps='sphere')

    y = icedata[icedata.dims[1]]
    x = icedata[icedata.dims[2]]

    dx = np.empty((icedata.shape[1], icedata.shape[2]))*np.nan
    dy = np.empty((icedata.shape[1], icedata.shape[2]))*np.nan

    for i in range(len(x)-1):
        for j in range(len(y)):
            _,_,dx[j,i] = g.inv(x[i].values, y[j].values, x[i+1].values, y[j].values)
    for j in range(len(y)):
        _,_,dx[j, -1] = g.inv(x[-1].values, y[j].values,x[0].values, y[j].values)
    
    for i in range(len(x)):
        for j in range(len(y)-1):
            _,_,dy[j,i] = g.inv(x[i].values, y[j].values, x[i].values, y[j+1].values)
    for i in range(len(x)):
        dy[-1, i] = dy[-2, i]
    
    areadata = xr.DataArray(
        data = dx*dy,
        dims=icedata.dims[1:],
        coords={list(icedata.coords)[1]: icedata[list(icedata.coords)[1]], list(icedata.coords)[2]:icedata[list(icedata.coords)[2]]}
    )
    return areadata

def match_unmatching_grid(icedata, areadata):
    newareadata = xr.DataArray(
        data=np.empty(icedata.shape[1:])*np.nan,
        dims=icedata.dims[1:],
        coords={list(icedata.coords)[1]: icedata[list(icedata.coords)[1]], list(icedata.coords)[2]:icedata[list(icedata.coords)[2]]}
    )
    return newareadata + areadata

def set_nan_to_zero(icedata):
    ## for E3SM-2-0
    d0 = icedata.fillna(0)
    mlat = icedata.idxmin(icedata.dims[1])  ## find the ice extent edge
    a = d0.where(d0[icedata.dims[1]]>=mlat).where(d0 == 0)
    newice = xr.where((a>=0)|(icedata>=0), d0, np.nan)
    return newice

def set_zero_to_nan(icedata):
    ## for GISS, INM-CM4-8
    dnan = icedata.where(icedata>0)
    mlat = dnan.idxmin(dnan.dims[len(dnan.dims)-2])  ## find the ice extent edge
    a = icedata.where(icedata[icedata.dims[len(icedata.dims)-2]]>=mlat).where(icedata == 0)
    newice = xr.where((a>=0)|(dnan>=0), icedata, np.nan)
    return newice

def modify_area_grid_to_ice_grid(icedata, areadata):
    ## only for 'CAS-ESM2-0'
    ## area lon 0~359 ; ice lon 1~360 
    ## first modify the lon in areadata
    a = areadata.isel({areadata.dims[1]:0})  # select lon=0
    a['lon'] = areadata[areadata.dims[1]][-1].values+1  # resign lon = 360 to lon=0 
    b = areadata.sel({areadata.dims[1]:slice(1, None)})  # select lon = 1~359
    new_area = xr.concat([b, a], dim=areadata.dims[1])  # combine 1~359 & 360
    area = copy_coords_xy(icedata, new_area, changename=True)
    return area


def flip_y(ds):
    ## for MPI 
    new_y = np.flip(ds[ds.dims[len(ds.dims)-2]])
    ds = ds.reindex({ds.dims[len(ds.dims)-2]: new_y})
    dsnew = ds.assign_coords(
        {
            ds.dims[len(ds.dims)-2]: range(0, len(ds[ds.dims[len(ds.dims)-2]]))
        }
    )
    return dsnew
