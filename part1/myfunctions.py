import xarray as xr
import pandas as pd
import numpy as np

import os
import pickle
import gcsfs
import glob

from pyproj import Geod

def ispickleexists(n, p0):
    p = p0 + n + '.pickle'
    if os.path.exists(p):
        # print('    [o] {} exists.'.format(p))
        return True
    else:
        return False

def openpickle(n, p0):
    p = p0 + n + '.pickle'
    d = pd.read_pickle(p)
    # with open(p, 'rb') as df:
    #     d = pickle.load(df)
    return d

def savepickle(n, p0, sf):
    p = p0 + n + '.pickle'
    with open(p, 'wb') as wf:
        pickle.dump(sf, wf, pickle.HIGHEST_PROTOCOL)

def open_from_cloud(link):
    gcs = gcsfs.GCSFileSystem(token='anon')
    mapper = gcs.get_mapper(link)
    ## open it using xarray and zarr
    ds = xr.open_zarr(mapper, consolidated=True)
    if 'type' in ds.coords:
        ds = ds.reset_coords('type', drop = True)
    return ds

def open_nc(mf):
    if len(mf)>50:
        ds = xr.open_mfdataset(mf[0], use_cftime=True)
        for i in range(1, len(mf)):
            ds0 = xr.open_mfdataset(mf[i], use_cftime=True)
            ds = xr.concat([ds, ds0], dim="time")
    else:
        ds = xr.open_mfdataset(mf, use_cftime=True)
    if 'type' in ds.coords:
        ds = ds.reset_coords('type', drop = True)
    return ds

def select_month(da, n):
    return da.isel(time=(da.time.dt.month == n))

def open_nc_month(mf, month):
    if len(mf)>50:
        ds = xr.open_mfdataset(mf[0], use_cftime=True)
        ds = select_month(ds, month)
        for i in range(1, len(mf)):
            ds0 = xr.open_mfdataset(mf[i], use_cftime=True)
            ds0 = select_month(ds0, month)
            ds = xr.concat([ds, ds0], dim="time")
        # ds = xr.open_mfdataset(mf, use_cftime=True)
    else:
        ds = xr.open_mfdataset(mf, use_cftime=True)
        ds = select_month(ds, month)
    if 'type' in ds.coords:
        ds = ds.reset_coords('type', drop = True)
    return ds

def read_siconc(p_nc, name, selected_month):
    data_path = p_nc + 'siconc_' + '*' + name + '_piControl_' + '*' + '.nc'
    matching_files = glob.glob(data_path)
    if len(matching_files) > 0:
        if selected_month:
            ds = open_nc_month(matching_files, selected_month) 
        else:
            ds = open_nc(matching_files)
    else:
        raise ValueError("    [x] no ice data.")
    return ds    

def get_latlon(datapd, i, ds0, newlatlon=False, nolatlon=False):
    if not nolatlon:
        nolatlon = pd.isna(datapd.at[i, 'latname'])
    if nolatlon:
        if newlatlon:
            latname, lonname = newlatlon
        else:
            latname = datapd.at[i, 'yname']
            lonname = datapd.at[i, 'xname']
        x = ds0[lonname]
        y = ds0[latname]
        newlon, newlat = np.meshgrid(x, y)
        dlon = xr.DataArray(newlon, dims={datapd.at[i, 'yname']:y.values, datapd.at[i, 'xname']:x.values})
        dlat = xr.DataArray(newlat, dims={datapd.at[i, 'yname']:y.values, datapd.at[i, 'xname']:x.values})
    else:
        if newlatlon:
            latname, lonname = newlatlon
        else:
            latname = datapd.at[i, 'latname']
            lonname = datapd.at[i, 'lonname']
        dlat = ds0[latname].load()
        dlon = ds0[lonname].load()
        dlat = dlat.where(dlat < 91).where(dlat>-91)
        dlon = dlon.where(dlon < 361).where(dlon>-361)  
    if 'time' in dlat.coords:
        dlat = dlat.reset_coords('time', drop = True)
        dlon = dlon.reset_coords('time', drop = True)
    if 'time_bounds' in dlat.coords:
        dlat = dlat.reset_coords('time_bounds', drop = True)
        dlon = dlon.reset_coords('time_bounds', drop = True)
    return dlat, dlon

def get_south(ds, datapd, i, southlat = -40):
    if time in ds.dims:
        dlat, lon = get_latlon(datapd, i, ds.isel(time=0))
    else:
        dlat, lon = get_latlon(datapd, i, ds)
    da_south = ds.where(dlat <= southlat, drop=True)
    return da_south

def read_areacello(p_nc, name, var):
    # Read cell area nc file
    grid_path = p_nc + var + '_Ofx_' + name + '_*' + '.nc'
    matching_gfiles = glob.glob(grid_path)
    if len(matching_gfiles)==1:
        dsg = xr.open_mfdataset(matching_gfiles)
    else:
        raise ValueError("    [x] cell area data error.")
    return dsg

def create_new_ds(ice, area, lat, lon):
    area_data = area.values
    try:
        new_ds = xr.Dataset(
            data_vars = {
                'siconc': ice,
                'areacello': (lat.dims, area_data),
                'newlat': lat,
                'newlon': lon,
            }
        )
        return new_ds
    except Exception as error:
        print("    An exception occurred:", error)

def calculate_area_xy(ds):
    g = Geod(ellps='sphere')
    icedata = ds.siconc
    ybnds_name = icedata.dims[1] + '_bnds'
    xbnds_name = icedata.dims[2] + '_bnds'
    if ybnds_name in ds.data_vars:
        ybnds = ds[ybnds_name].values
        xbnds = ds[xbnds_name].values
    elif ybnds_name in ds.coords:
        ybnds = ds[ybnds_name].values
        xbnds = ds[xbnds_name].values
    else:
        print(ds.data_vars)
        raise Error("    no coords bnds")

    y = ds[icedata.dims[1]].values
    x = ds[icedata.dims[2]].values

    dx = np.empty((icedata.shape[1], icedata.shape[2]))*np.nan
    dy = np.empty((icedata.shape[1], icedata.shape[2]))*np.nan

    for i in range(len(x)):
        for j in range(len(y)):
            _,_,dx[j,i] = g.inv(xbnds[i,0], y[j], xbnds[i,1], y[j])          
    
    for j in range(len(y)):
        _,_,dy[j,:] = g.inv(x[0], ybnds[j, 0], x[0], ybnds[j, 1])
    
    areadata = xr.DataArray(
        data = dx*dy,
        dims = icedata.dims[1:],
        coords={list(icedata.dims)[-1]: icedata[list(icedata.dims)[-1]], list(icedata.dims)[-2]:icedata[list(icedata.dims)[-2]]}
    )
    return areadata

# def calculate_area_latlon(lat, lon):
    

def newxy_fmissingxy(dx, dy):
    dxv = dx.where((dx>-361) & (dx<361)).values
    dyv = dy.where((dy>-91) & (dy<91)).values
    newx0 = dxv[~np.isnan(dxv).any(axis=1)][0]
    newy0 = dyv[:, ~np.isnan(dyv).any(axis=0)][:,0]
    newx, newy = np.meshgrid(newx0, newy0)
    x = np.where(np.isnan(dx), newx, dxv)
    y = np.where(np.isnan(dy), newy, dyv)
    dsx = xr.DataArray(x, dims = dx.dims, coords = dx.coords)
    dsy = xr.DataArray(x, dims = dx.dims, coords = dx.coords)
    return dsx, dsy

# def get_new_xy(d, datapd, j, newlonlat = 0):
#     if newlonlat:
#         newlon = d[newlonlat[0]]
#         newlat = d[newlonlat[1]]
#     else:
#         if pd.isna(datapd.at[j, 'latname']):
#             pltx0 = d[datapd.at[j, 'xname']]
#             plty0 = d[datapd.at[j, 'yname']]
#             newlon0, newlat0 = np.meshgrid(pltx0, plty0)
#             newlon = xr.DataArray(newlon0, dims={datapd.at[j, 'yname']:plty0.values, datapd.at[j, 'xname']:pltx0.values})
#             newlat = xr.DataArray(newlat0, dims={datapd.at[j, 'yname']:plty0.values, datapd.at[j, 'xname']:pltx0.values})
#         else:
#             newlon = d[datapd.at[j, 'lonname']]
#             newlat = d[datapd.at[j, 'latname']]
#     if len(np.shape(newlat)) > 2:
#         newlon = newlon.isel(time = 0)
#         newlat = newlat.isel(time = 0)
#     return newlon, newlat



# def add_newlatlon(datapd, i, ds):
#     dlat, lon = get_latlon(datapd, i, ds.isel(time=0))
#     newlon, newlat = newxy_fmissingxy(dlon, dlat)
#     ds['newlat'] = (dlat.dims, newlat)
#     ds['newlon'] = (dlon.dims, newlon)
#     return ds






# # def get_sepsouth(mf, datapd, i, southlat = -40):
# #     ds = xr.open_mfdataset(mf, use_cftime=True)
# #     da = ds.siconc
# #     if 'type' in da.coords:
# #         da = da.reset_coords('type', drop = True)
# #     da_sep = da.isel(time=(da.time.dt.month == 9))
# #     da_south = get_south(da_sep, datapd, i, southlat)
# #     return da_south

# # def openicenc(p0, datapd, i, southlat = -40):
# #     data_path = p0 + datapd.at[i, 'source_id'] + '_piControl_' + datapd.at[i, 'member_id'] + '*' + datapd.at[i, 'grid_label'] + '*' + '.nc'
# #     if datapd.at[i, 'source_id']  == 'NorESM2-LM' or datapd.at[i, 'source_id'] == 'NorESM2-MM':
# #         data_path = p0 + datapd.at[i, 'source_id'] + '_piControl_' + datapd.at[i, 'member_id'] + '*' + '.nc'
# #     matching_files = glob.glob(data_path)
# #     if len(matching_files)>200:    
# #         for t in range(len(matching_files)):
# #             da_south = get_sepsouth(matching_files[t], datapd, i, southlat = southlat)
# #             if t == 0:
# #                 da_save = da_south.load()
# #             else:
# #                 da_save0 = da_south.load()
# #                 da_save = xr.concat([da_save, da_save0], dim="time")
# #         da_s = da_save
# #     else:
# #         da_south = get_sepsouth(matching_files, datapd, i, southlat = southlat)
# #         da_s = da_south.load()
# #     return da_s

def drop_coords(ds):
    for v in ds.coords:
        if v not in ds.dims:
            ds = ds.drop_vars(v)
    return ds
    
def copy_coords(copyfrom, copyto, latname, lonname):
    newd = copyto.assign_coords(
        {
            latname : copyfrom[latname],
            lonname : copyto[lonname],
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

def copy_xy(data_copyfrom, data_copyto):
    data_copyto = rename_xy(data_copyfrom, data_copyto)
    data_copyto[data_copyto.dims[len(data_copyto.dims)-1]] = data_copyfrom[data_copyfrom.dims[len(data_copyfrom.dims)-1]].values
    data_copyto[data_copyto.dims[len(data_copyto.dims)-2]] = data_copyfrom[data_copyfrom.dims[len(data_copyfrom.dims)-2]].values
    return data_copyto


# def replace_missingxg(d, datapd, j, newlonlat = 0):
#     if newlonlat:
#         lon = d[newlonlat[0]]
#         lat = d[newlonlat[1]]
#     else:
#         lon = d[datapd.at[j, 'lonname']]
#         lat = d[datapd.at[j, 'latname']]
#     if 'time' in lat.dims:
#         lon = lon.isel(time = 0)
#         lat = lat.isel(time = 0)
#     newlon, newlat = newxy_fmissingxy(lon, lat)
#     newd = d.assign_coords(
#         {
#             datapd.at[j, 'lonname']: (lon.dims, newlon),
#             datapd.at[j, 'latname']: (lat.dims, newlat),
#         }
#     )
#     return newd



# def match_unmatching_grid(icedata, areadata):
#     newareadata = xr.DataArray(
#         data=np.empty(icedata.shape[1:])*np.nan,
#         dims=icedata.dims[1:],
#         coords={list(icedata.coords)[1]: icedata[list(icedata.coords)[1]], list(icedata.coords)[2]:icedata[list(icedata.coords)[2]]}
#     )
#     return newareadata + areadata

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

# def modify_area_grid_to_ice_grid(icedata, areadata):
#     ## only for 'CAS-ESM2-0'
#     ## area lon 0~359 ; ice lon 1~360 
#     ## first modify the lon in areadata
#     a = areadata.isel({areadata.dims[1]:0})  # select lon=0
#     a['lon'] = areadata[areadata.dims[1]][-1].values+1  # resign lon = 360 to lon=0 
#     b = areadata.sel({areadata.dims[1]:slice(1, None)})  # select lon = 1~359
#     new_area = xr.concat([b, a], dim=areadata.dims[1])  # combine 1~359 & 360
#     area = copy_coords_xy(icedata, new_area, changename=True)
#     return area


# def flip_y(ds):
#     ## for MPI 
#     new_y = np.flip(ds[ds.dims[len(ds.dims)-2]])
#     ds = ds.reindex({ds.dims[len(ds.dims)-2]: new_y})
#     dsnew = ds.assign_coords(
#         {
#             ds.dims[len(ds.dims)-2]: range(0, len(ds[ds.dims[len(ds.dims)-2]]))
#         }
#     )
#     return dsnew

