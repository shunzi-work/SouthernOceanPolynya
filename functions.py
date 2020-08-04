# -*- coding: UTF-8 -*-
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
        naninds = np.where(np.isnan(dens_diff))[0]
        if len(naninds) > 0:
            nanindex = naninds[0]
        else:
            nanindex = len(depths)
        dens_diff_drop = dens_diff[:nanindex]
        if np.all(np.diff(dens_diff_drop) > 0):
            if dens_diff_drop[-1] > 0: 
                mld = np.interp(0, dens_diff, depths)
            else:
                mld = depths[nanindex-1]
        else: # A tolerance threshold, dens < dens(10m)-0.01 or dens > dens(10m) + 0.03
            new_diff = dens_diff[1:]
            nthr_index = np.where(np.abs(new_diff + 0.02) > 0.02)[0]
            if len(nthr_index) == 0:
                mld = depths[nanindex-1]
            else:
                nind = nthr_index[0] + 2
                dens_diff_inthr = dens_diff[:nind]
                if np.all(dens_diff_inthr<0):
                    mld = np.nan
                else:
                    mld = np.interp(0, dens_diff_inthr, depths[:nind])
    return mld
