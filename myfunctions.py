import xarray as xr
import pandas as pd
import numpy as np

import os
import pickle
import gcsfs
import glob

from pyproj import Geod
from scipy import ndimage
from skimage.segmentation import flood_fill

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
    if (ybnds_name in ds.data_vars) or (ybnds_name in ds.coords):
        ybnds = ds[ybnds_name].values
        xbnds = ds[xbnds_name].values
    else:
        raise TypeError("No x/y bnds")

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


def newxy_fmissingxy(dx, dy):
    dxv = dx.where((dx>-361) & (dx<361)).values
    dyv = dy.where((dy>-91) & (dy<91)).values
    newx0 = dxv[~np.isnan(dxv).any(axis=1)][0]
    newy0 = dyv[:, ~np.isnan(dyv).any(axis=0)][:,0]
    newx, newy = np.meshgrid(newx0, newy0)
    x = np.where(np.isnan(dx), newx, dxv)
    y = np.where(np.isnan(dy), newy, dyv)
    dsx = xr.DataArray(x, dims = dx.dims, coords = dx.coords)
    dsy = xr.DataArray(y, dims = dx.dims, coords = dx.coords)
    return dsx, dsy


def find_first_non_nan_row(da):
    """
    Finds the first row in each column where the value is not NaN.
    Returns:
        An xarray DataArray containing the coordinates of the first non-NaN value in each column.
    """
    # Find non-NaN values
    da = da.where(da > 0)
    not_nan = ~da.isnull()
    # Get the indices of the first non-NaN value 
    first_non_nan_indices = not_nan.argmax(da.dims[-2])
    # Get the corresponding coordinate values
    first_non_nan_rows = da[da.dims[-2]][first_non_nan_indices]
    return first_non_nan_rows


def detect_polynya(daice, daarea, ice_threshold, area_threshold, buffering = 0):
    s = ndimage.generate_binary_structure(2,2)
    da_masked = xr.DataArray(np.nan*np.empty_like(daice), dims = daice.dims, coords = daice.coords)
    for year in daice.time:
        ice0 = daice.sel(time = year)
        icenew = ice0 <= ice_threshold
        ice = xr.where(np.isnan(ice0), True, icenew)   # get rid of "coastal polynya" 
        if buffering: # buffering (more layers)
            b = 0
            while b<buffering:
                ice = xr.where(ice0[ice0.dims[-2]] == find_first_non_nan_row(ice0), True, ice) 
                b+=1
        labeled_image, num_features = ndimage.label(ice, structure = s)
        if num_features < 2:
            continue
        mask = np.zeros_like(labeled_image)
        for i in range(1, num_features+1):
            area = daarea.where(labeled_image == i).sum()/1e6  # m2 -> km2
            if (area > area_threshold[0]) and (area < area_threshold[1]):  # the area of open 'polynya' within the sea ice extent is small
                mask[labeled_image == i] = 1
        da_masked.loc[year] = xr.where(mask, ice0, np.nan)
    return da_masked #.mean('time'), da_masked.count('time')


def count_polynya_area(ds, ice_threshold, area_threshold, b=0, re=False):
    masked = detect_polynya(ds.siconc, ds.areacello, ice_threshold, area_threshold, b)
    polynya_count = masked.count('time')
    if re:
        polynya_count = polynya_count.where(polynya_count >= re)
    return ds.areacello.where(polynya_count > 0).sum().values.item()


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
            data_copyto.dims[-1]:data_copyfrom.dims[-1], 
            data_copyto.dims[-2]:data_copyfrom.dims[-2], 
        }
    )
    return new


def copy_xy(data_copyfrom, data_copyto):
    data_copyto = rename_xy(data_copyfrom, data_copyto)
    data_copyto[data_copyto.dims[-1]] = data_copyfrom[data_copyfrom.dims[-1]].values
    data_copyto[data_copyto.dims[-2]] = data_copyfrom[data_copyfrom.dims[-2]].values
    return data_copyto


def set_land_to_nan(ds):
    ## for GISS, INM-CM4-8
    ice = ds.siconc
    for t in ice.time:
        ice0 = ice.sel(time = t)
        dffill = flood_fill(ice0.values, (0, 0), np.nan, tolerance=0)
        if ~np.isnan(dffill[-1,-1]):
            break
    ice_mask = xr.DataArray(
        data = dffill,
        dims = ice0.dims, 
        coords=ice0.coords
    )
    newice = ice.where(~ice_mask.isnull())
    ds['siconc'] = newice
    return ds

def set_ocean_to_zero(ds):
    ## for E3SM-2-0
    dsnew = ds.siconc.fillna(0)
    ds['siconc'] = dsnew
    newice = set_land_to_nan(ds)
    return newice

def shift_x(da):
    ## only for 'CAS-ESM2-0'
    ## area lon 0~359 ; ice lon 1~360 
    ## first modify the lon in areadata
    ## NOTE: no need, those two actually match
    a = da.isel({da.dims[1]:0})  # select lon=0
    a[da.dims[1]] = da[da.dims[1]][0].values+360  # resign lon = 360 to lon=0 
    b = da.sel({da.dims[1]:slice(1, None)})  # select lon = 1~359
    new = xr.concat([b, a], dim=da.dims[1])  # combine 1~359 & 360
    return new

def change_start_x(ds, newx):
    lon360 = (ds.newlon+360)%360 - newx
    if np.abs(lon360[0,0])>5:
        diffmin = np.argmin(np.abs(lon360[0,:]).values)
        part_a = ds.isel({lon360.dims[-1]: slice(diffmin, None)})
        part_b = ds.isel({lon360.dims[-1]: slice(0, diffmin)})
        newds = xr.concat([part_a, part_b], dim = lon360.dims[-1])
        return newds
    else:
        return ds

def flip_y(ds):
    dsa = ds.reindex({ds.dims[-2]: ds[ds.dims[-2]][::-1]})
    dsnew = dsa.assign_coords({dsa.dims[-2]: ds[ds.dims[-2]]})
    return dsnew