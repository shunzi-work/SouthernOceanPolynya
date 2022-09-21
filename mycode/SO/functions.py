# Functions for CMIP6 model output data (xarray data) Processing

import numpy as np
import xesmf as xe
import xarray as xr

def smow(t):
    """
    Modified from python-seawater
    Density of Standard Mean Ocean Water (Pure Water) using EOS 1980.
    Parameters
    ----------
    t : array_like
        temperature [℃ (ITS-90)]
    Returns
    -------
    dens(t) : array_like
              density  [kg m :sup:`3`]
    Examples
    --------
    >>> # Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22.
    >>> import seawater as sw
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> sw.smow(t)
    array([ 999.842594  ,  999.842594  ,  995.65113374,  995.65113374,
            999.842594  ,  999.842594  ,  995.65113374,  995.65113374])
    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf
    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
       of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
       doi:10.1016/0198-0149(81)90122-9
    """
    a = (999.842594, 6.793952e-2, -9.095290e-3, 1.001685e-4, -1.120083e-6,
         6.536332e-9)
    T68 = t * 1.00024
    return (a[0] + (a[1] + (a[2] + (a[3] + (a[4] + a[5] * T68) * T68) * T68) *
            T68) * T68)

def dens0(s, t):
    """
    Modifed from python-seawater
    Density of Sea Water at atmospheric pressure.
    Parameters
    ----------
    s(p=0) : array_like
             salinity [psu (PSS-78)]
    t(p=0) : array_like
             temperature [℃ (ITS-90)]
    Returns
    -------
    dens0(s, t) : array_like
                  density  [kg m :sup:`3`] of salt water with properties
                  (s, t, p=0) 0 db gauge pressure
    Examples
    --------
    >>> # Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22
    >>> import seawater as sw
    >>> from seawater.library import T90conv
    >>> s = [0, 0, 0, 0, 35, 35, 35, 35]
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> sw.dens0(s, t)
    array([  999.842594  ,   999.842594  ,   995.65113374,   995.65113374,
            1028.10633141,  1028.10633141,  1021.72863949,  1021.72863949])
    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
       computation of fundamental properties of seawater. UNESCO Tech. Pap. in
       Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
       http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf
    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere
       equation of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
       doi:10.1016/0198-0149(81)90122-9
    """
    T68 = t * 1.00024
    b = (8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9)
    c = (-5.72466e-3, 1.0227e-4, -1.6546e-6)
    d = 4.8314e-4
    return (smow(t) + (b[0] + (b[1] + (b[2] + (b[3] + b[4] * T68) * T68) *
            T68) * T68) * s + (c[0] + (c[1] + c[2] * T68) * T68) * s *
            s ** 0.5 + d * s ** 2)

def func_mld(dens_diff, depths):
    '''
    Calculating the mixed layer depth based on the constant potential density 
    difference criterion.
    MLD = depth where(sigma[mld] = sigma[10] + 0.03 kg/m3). 
    (0.03 kg/m3 ~= 0.03 psu)
    ----------
    Parameters
    dens_diff: Data array of density difference [density - density(at 10m) - 0.03]
    depths:    Data array of depth 
    ----------
    References
    .. [1] de Boyer Montégut, C., G. Madec, A. S. Fischer, A. Lazar, and 
    D. Iudicone, 2004: Mixed layer depth over the global ocean: an examination
    of profile data and a profile-based climatology. J. Geophys. Res., 109, 
    C12003. doi:10.1029/2004JC002378
    '''
    if np.isnan(dens_diff[0]):
        mld = np.nan
    elif dens_diff[0] >= 0:
        mld = np.nan
    else:
        nthr_index = np.where(dens_diff > 0)[0]
        if len(nthr_index) == 0:
            naninds = np.where(np.isnan(dens_diff))[0]
            if len(naninds) > 0:
                nanindex = naninds[0]
            else:
                nanindex = len(depths)
            mld = depths[nanindex-1]
        else:
            nind = nthr_index[0] + 1
            mld = np.interp(0, dens_diff[:nind], depths[:nind])
    return mld

def func_regrid(ds, ds_out, gr_method = 'bilinear', reuse = False, clear = True):
    dsr = xe.Regridder(ds, ds_out, gr_method, periodic = True, reuse_weights = reuse, ignore_degenerate = True)
    dsr._grid_in = None
    dsr._grid_out = None
    dsr_out = dsr(ds)
    if clear:
        dsr.clean_weight_file()
    return dsr_out

def sel_time(ds, start_year, end_year, month = None):
    ds = ds.isel(time = slice((start_year-1)*12, end_year*12))
    if month:
        ds = list(ds.groupby("time.month"))[month-1][-1]
    return ds

def xr_func_mld(dens):
    dens10 = dens.interp(lev = 10, method = 'linear')  # density at 10m
    dens_diff = dens - dens10 - 0.03               # density differences 
    mld = xr.apply_ufunc(
        func_mld, 
        dens_diff.chunk({"time":25,"lon":45,"lat":45}),  
        dens_diff.lev, 
        input_core_dims = [["lev"], ["lev"]], 
        vectorize = True,
        dask = "parallelized",
        output_dtypes = [dens_diff.lev.dtype])
    return mld

def select_conv_area_data(ds, conv_area):
    ds = ds.copy()
    ds = ds.where(conv_area > 0)
    ds = ds.groupby('time.year').mean(dim = 'time', skipna = True)
    ds = ds.mean(dim = 'lon', skipna = True).mean(dim = 'lat', skipna = True)
    return ds

def select_region(ds, region, reg_dict):
    return ds.sel(lon = reg_dict[region]['lon'], lat = reg_dict[region]['lat'], method = "nearest")

def add_region_attrs(ds, region, reg_dict):
    ds.attrs['region_name'] = reg_dict[region]['name']
    ds.attrs['conv_index'] = reg_dict[region]['convind']
    ds.attrs['rep_lon'] = reg_dict[region]['rlon']
    ds.attrs['rep_lat'] = reg_dict[region]['rlat']
    return ds

def func_conv_data(datasets, region, start_year, end_year, reg_dict, month_no = 9, regrid = None):
    da_t = sel_time(datasets['thetao'].thetao, start_year, end_year)
    da_s = sel_time(datasets['so'].so, start_year, end_year)
    da_dens = sel_time(dens0(da_s, da_t), start_year, end_year, month = month_no)
    if not regrid == None:
        da_dens = func_regrid(da_dens, regrid)
        da_t = func_regrid(da_t, regrid)
        da_s = func_regrid(da_s, regrid)
    da_dens = select_region(da_dens, region, reg_dict)
    da_t = select_region(da_t, region, reg_dict)
    da_s = select_region(da_s, region, reg_dict)
    da_mld = xr_func_mld(da_dens)
    conv = xr.where(da_mld >= reg_dict[region]['convdepth'], da_mld, np.nan)
    conv_area = conv.mean(dim = 'time', skipna = True)
    da_t_conv = select_conv_area_data(da_t, conv_area)
    da_s_conv = select_conv_area_data(da_s, conv_area)
    conv_area = add_region_attrs(conv_area, region, reg_dict)
    da_t_conv = add_region_attrs(da_t_conv, region, reg_dict)
    da_s_conv = add_region_attrs(da_s_conv, region, reg_dict)
    ind_t = da_t_conv.interp(lev = da_t_conv.conv_index)
    conv_ind = (ind_t - ind_t.mean('year'))/ind_t.std('year') * -1
    return conv_area, da_t_conv, da_s_conv, conv_ind


def func_select_oceanindex(lat, lon, ds):
    '''
    Modified from https://stackoverflow.com/questions/58758480/xarray-select-nearest-lat-lon-with-multi-dimension-coordinates
    '''
    abslat = np.abs(ds.TLAT-lat)
    abslon = np.abs(ds.TLONG-lon)
    c = np.maximum(abslat, abslon) 
    latloc, lonloc = np.where(c == np.min(c))
    ocean_index = ds.REGION_MASK.sel(nlat=latloc[0], nlon=lonloc[0]).values - 0
    return ocean_index


def func_add_mask(ds):
    mask_ocean = 2 * np.ones((ds.dims['lev'], ds.dims['y'], ds.dims['x'])) * np.isfinite(ds.uo.isel(time=0))  
    mask_land = 1 * np.ones((ds.dims['lev'], ds.dims['y'], ds.dims['x'])) * np.isnan(ds.uo.isel(time=0))  
    mask_array = mask_ocean + mask_land
    ds.coords['mask'] = (('lev', 'y', 'x'), mask_array)
    return ds

def calc_dx_dy(ds, bnd = False):
    ''' This definition calculates the distance 
        between grid points that are in
        a latitude/longitude format.
        
        Using pyproj GEOD; different Earth Shapes 
        https://jswhit.github.io/pyproj/pyproj.Geod-class.html
        Common shapes: 'sphere', 'WGS84', 'GRS80'
        
        Accepts, 1D arrays for latitude and longitude
        
        Returns: dx, dy; 2D arrays of distances 
                       between grid points in the x and y direction in meters
        ------------
        modified from:
        https://github.com/Unidata/MetPy/issues/288#issuecomment-279481555
    '''
    from pyproj import Geod
    
    g = Geod(ellps='sphere')
    
    lon = ds.lon.values
    lat = ds.lat.values
    
    dx = np.empty(ds.lon.shape)
    dy = np.zeros(ds.lat.shape)
    
    for i in range(dx.shape[0]-1):
        for j in range(dx.shape[1]):
            _, _, dy[i,j] = g.inv(lon[i,j], lat[i,j], lon[i+1,j], lat[i+1,j])
    dy[i+1,:] = dy[i,:]
            
    for i in range(dx.shape[0]):
        for j in range(dy.shape[1]-1):
            _, _, dx[i,j] = g.inv(lon[i,j], lat[i,j], lon[i,j+1], lat[i,j+1])
    dx[:,j+1] = dx[:,j]
    
    return dx, dy