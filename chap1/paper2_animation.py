import numpy as np
import xesmf as xe
import xarray as xr
import pandas as pd

import gsw
import copy
import os
import pickle
import time

from xgcm import Grid
from pyproj import Geod
from scipy import stats

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors
import matplotlib.path as mpath
import matplotlib.animation as animation
# import cmocean

from matplotlib.ticker import ScalarFormatter

def load_data(dataname):
    # load data either from save pickle file 
    data_path = 'newds/'+ dataname + '.pickle'
    if os.path.exists(data_path):
        with open(data_path, 'rb') as f:
            out = pickle.load(f)
    return out

conv_index_rs = load_data('conv_index_rs')
# conv_index_ws = load_data('conv_index_ws')
# conv_index_ws_400 = load_data('conv_index_ws_400')
Psi_r_ann = load_data('Psi_r_ann')
# Psi_w_ann = load_data('Psi_w_ann')

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

fig = plt.figure(figsize=(10, 8)) 

dz = load_data('dz')

pltdata1 = load_data('newds_sep_south_400_1000_so')/dz.sel(lev=slice(401,1000)).sum('lev')
pltdata1 = pltdata1 - pltdata1.mean('time')
pltdata2 = load_data('newds_sep_south_400_1000_thetao')/dz.sel(lev=slice(401,1000)).sum('lev')
pltdata2 = pltdata2 - pltdata2.mean('time')

def plot_each_year(fig, i):
    plt.clf()
    plt.suptitle('400 m - 1000 m')

    axl = fig.add_subplot(2,1,1)
    axl.plot(conv_index_rs, label = 'Ross index')
    Psi_plot = (Psi_r_ann - Psi_r_ann.mean())/Psi_r_ann.std()
    axl.plot(Psi_plot, label = 'Ross gyre')
    axl.legend(ncol=2, loc='lower left')
    axl.axvline(x=i+100, c = 'r')

    axr = fig.add_subplot(2,2,3, projection=ccrs.SouthPolarStereo())
    axr.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())
    im = axr.pcolormesh(pltdata1.lon, pltdata1.lat, pltdata1.isel(time=i), 
                    transform=ccrs.PlateCarree(), vmin=-1, vmax=1, cmap=plt.cm.RdBu_r)
    axr.add_feature(cfeature.LAND, zorder=1)#, color='0.8')
    axr.add_feature(cfeature.COASTLINE, linewidth=0.2)
    axr.set_boundary(circle, transform=axr.transAxes)
    cb = fig.colorbar(im)
    axr.set_title('salt') 
    ax3 = fig.add_subplot(2,2,4, projection=ccrs.SouthPolarStereo())

    ax3.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())
    im = ax3.pcolormesh(pltdata2.lon, pltdata2.lat, pltdata2.isel(time=i), 
                    transform=ccrs.PlateCarree(), vmin=-3, vmax=3, cmap=plt.cm.RdBu_r)
    ax3.add_feature(cfeature.LAND, zorder=1)#, color='0.8')
    ax3.add_feature(cfeature.COASTLINE, linewidth=0.2)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    cb = fig.colorbar(im)
    ax3.set_title('temperature')
    return fig


for i in range(len(pltdata1.time)):
    fig = plot_each_year(fig, i)
    savename = 'paper2fig/img' + str(i+1) + '.png'
    fig.savefig(savename)
    print(i)
    # plt.show()




'''
## code copied from pangeo

# !pip install cmocean

import numpy as np
import xesmf as xe
import xarray as xr
import pandas as pd

import gsw
import copy
import os
import pickle
import time

from xgcm import Grid
from pyproj import Geod
from scipy import stats

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as colors
import matplotlib.path as mpath
# import cmocean

from matplotlib.ticker import ScalarFormatter

### fetch data from google cloud
from intake import open_esm_datastore
col = open_esm_datastore("https://storage.googleapis.com/cmip6/pangeo-cmip6.json")

## salt and temp
dslist = col.search(variable_id = ['so', 'thetao'], table_id = 'Omon', 
                    experiment_id = 'piControl', source_id = 'GFDL-CM4', grid_label = 'gn')
ds = dslist['CMIP.NOAA-GFDL.GFDL-CM4.piControl.Omon.gn'].to_dask()
ds = ds.squeeze('member_id').reset_coords('member_id', drop = True).squeeze('dcpp_init_year').reset_coords('dcpp_init_year', drop = True)

import zarr
import gcsfs
# this only needs to be created once
gcs = gcsfs.GCSFileSystem(token='anon')

# area and volume 
mapper_area = gcs.get_mapper('gs://cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/piControl/r1i1p1f1/Ofx/areacello/gn/v20180701')
mapper_vol = gcs.get_mapper('gs://cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/piControl/r1i1p1f1/Omon/volcello/gn/v20180701')

ds_area = xr.open_zarr(mapper_area, consolidated=True)
ds_vol = xr.open_zarr(mapper_vol, consolidated=True)


# wmo
mapper_wmo = gcs.get_mapper('gs://cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/piControl/r1i1p1f1/Omon/wmo/gn/v20180701/')
ds_wmo = xr.open_zarr(mapper_wmo, consolidated=True)


# u and v
mapper_uo = gcs.get_mapper('gs://cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/piControl/r1i1p1f1/Omon/uo/gn/v20180701/')
mapper_vo = gcs.get_mapper('gs://cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/piControl/r1i1p1f1/Omon/vo/gn/v20180701/')

ds_uo = xr.open_zarr(mapper_uo, consolidated=True)
ds_vo = xr.open_zarr(mapper_vo, consolidated=True)

# wind 
mapper_tauuo = gcs.get_mapper('gs://cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/piControl/r1i1p1f1/Omon/tauuo/gn/v20180701')
mapper_tauvo = gcs.get_mapper('gs://cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/piControl/r1i1p1f1/Omon/tauvo/gn/v20180701')

ds_tauuo = xr.open_zarr(mapper_tauuo, consolidated=True)
ds_tauvo = xr.open_zarr(mapper_tauvo, consolidated=True)


# ocean basin region
mapper_basin = gcs.get_mapper('gs://cmip6/CMIP6/CMIP/NOAA-GFDL/GFDL-CM4/piControl/r1i1p1f1/Ofx/basin/gn/v20180701')
ds_basin = xr.open_zarr(mapper_basin, consolidated=True)


# ugrid = xr.open_dataset('ocean_static_0p25_ugrid.nc').load()
vgrid = xr.open_dataset('ocean_static_0p25_vgrid.nc').load() 

def calc_dx_dy_bnd(ds):
    g = Geod(ellps='sphere')
    
    lonb = ds.lon_bnds.values
    latb = ds.lat_bnds.values
    
    dx = np.zeros((lonb.shape[0], lonb.shape[1]))
    dy = np.zeros((lonb.shape[0], lonb.shape[1]))
    
    for i in range(dx.shape[0]):
        for j in range(dx.shape[1]):
            _, _, dy[i,j] = g.inv(lonb[i,j,0], latb[i,j,0], lonb[i,j,3], latb[i,j,3])
            
    for i in range(dx.shape[0]):
        for j in range(dy.shape[1]):
            _, _, dx[i,j] = g.inv(lonb[i,j,0], latb[i,j,0], lonb[i,j,1], latb[i,j,1])
    
    return dx, dy

dx, __ = calc_dx_dy_bnd(ds)
dx = xr.DataArray(dx, coords=[vgrid.y, vgrid.x], dims=["y", "x"])
# dy = xr.DataArray(dy, coords=[U_soat_corrected.y, U_soat_corrected.x], dims=["y", "x"])


# uo_soat = ds_uo.uo.where((ds_basin.basin == 1) | (ds_basin.basin == 2))
vo_soat = ds_vo.vo.where((ds_basin.basin == 1) | (ds_basin.basin == 2))
vo_sopa = ds_vo.vo.where((ds_basin.basin == 1) | (ds_basin.basin == 3))
vo_sopi = ds_vo.vo.where((ds_basin.basin == 1) | (ds_basin.basin == 5) | (ds_basin.basin == 3))

dz = (ds.lev_bnds[:,1] - ds.lev_bnds[:,0]).load()

def cal_streamfunction(vo, dz, vgrid, dx):
    vdz = vo*dz
    V = vdz.sum('lev') - vdz.cumsum('lev')
    V_corrected = V.assign_coords(x = vgrid.x, y = vgrid.y)
    V_corrected = V_corrected.assign_coords(lon = vgrid.lon, lat = vgrid.lat)
    Vdx = V_corrected * dx /1e6
    Psi_max = Vdx.sum('x').max('lev')  #AABW
    Psi_min = Vdx.sum('x').min('lev')  #AAIW
    return Psi_max, Psi_min

Psi_soat_max, Psi_soat_min = cal_streamfunction(vo_soat, dz, vgrid, dx)
Psi_sopa_max, Psi_sopa_min = cal_streamfunction(vo_sopa, dz, vgrid, dx)
Psi_sopi_max, Psi_sopi_min = cal_streamfunction(vo_sopi, dz, vgrid, dx)


def load_data(dataname):
    # load data either from save pickle file or load data after computation
    data_path = dataname + '.pickle'
    if os.path.exists(data_path):
        with open(data_path, 'rb') as f:
            out = pickle.load(f)
    else:
        data = globals()[dataname]
        out = data.load()
        with open(data_path, 'wb') as f:
            pickle.dump(out, f, pickle.HIGHEST_PROTOCOL)
    return out


def MovingAverage(ds0, n):
    for i in range(len(ds0) - n):
        ds_avg = ds0[i:i+n].sum()/n
    return ds

def MovingAverage3D(ds):
    newds = []
    for yi in ds.y:
        ds0 = ds.sel(y=yi)
        ds_smooth = MovingAverage(ds0)
        newds.append(ds_smooth)
    
    
def lag_cor(x, y, lag):
    stat=[]
    x = (x-np.mean(x))/np.std(x)
    yind = (y-np.mean(y))/np.std(y)
    for i in range(2*lag+1):
        corr = np.corrcoef(yind[lag:len(x)-lag], x[i:len(x)-2*lag+i])
        stat.append(corr[1,0])
    return stat

def lag_cor_2D(x, y, lag):
    coors = np.zeros((len(x.lev), 2*lag+1))
    for l in range(len(x.lev)):
        x0 = x.isel(lev=l)
        stat = lag_cor(x0, y, lag)
        coors[l,:] = stat
    return coors



'''