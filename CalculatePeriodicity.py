###############
# 
# # Python script for calculating some properties 
# # reagarding polynya and convection each year
# 
###############

from myfunctions import *
import gc
from CalculateProperties import check_timerange, open_mld

import scipy.fft as fft
import scipy.stats as stats
import matplotlib.pyplot as plt
# filter some warning messages
import warnings
warnings.filterwarnings("ignore")

def calculate_alpha(ds):
    dsn = (ds - ds.mean())/ds.std()
    N = len(ds)
    acf = np.correlate(dsn, dsn, 'full')
    lag0 = acf[N-1]/N
    lag1 = acf[N]/(N-1)
    return lag1/lag0

def red_noise_spectrum(N, alpha):
    # contstruct expected red noise spectrum
    # use the Gilman et al. expression for the power spectrum of a red noise process
    rspec = np.empty(N)
    for h in np.arange(0, N):
        rspec[h] = (1.-alpha**2)/(1 -2 * alpha * np.cos(np.pi*h/(N-1))+alpha**2)
    return rspec

def cal_fft(ds):
    T = len(ds)
    freq = fft.rfftfreq(T, 1)
    periods = 1 / freq
    # compute power spectrum
    xf = fft.rfft(ds.values-ds.mean().item())
    pave = (2 / T * xf * xf.conj()).real
    # p = sig.welch(ds, window='hann', nperseg=T)
    # pave = p[1]
    # normalize the spectrum
    pave = pave/np.sum(pave)
    return periods, pave

def cal_sig_value(rspec, sig, dof=2, rdof=500):
    fstat = stats.f.ppf(sig,dof,rdof)
    return fstat*rspec

def plot_spectrum(ds, plot_title):
    alpha = calculate_alpha(ds)     #red noise lag-one autocorrelation
    periods, pave = cal_fft(ds)
    rspec = red_noise_spectrum(len(periods), alpha)
    
    # calculate significance using F-test
    spec99 = cal_sig_value(rspec, 0.99)
    spec90 = cal_sig_value(rspec, 0.90)
        
    # plot power spectrum and red noise spectra
    fig, ax = plt.subplots(figsize=(7,5))
    ax.set_xlabel('Period')
    ax.set_ylabel('Normalized Power')
    ax.set_title(plot_title)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter('{x:,.0f}')
    ax.plot(periods,pave,'-k', label = 'data')
    ax.plot(periods,rspec/np.sum(rspec),'-', label = 'red-noise fit', color = 'red')
    ax.plot(periods,spec99/np.sum(rspec),'--', label = '99% confidence', color = 'blue')
    ax.plot(periods,spec90/np.sum(rspec),'-.', label = '90% confidence', color = 'orange')
    ax.legend(bbox_to_anchor=(0., -0.02, 1., .102), loc='lower left', frameon=False,
                      ncols=4, mode="expand", borderaxespad=0.)
    fig.tight_layout()
    fig.savefig('fft_figures/' + plot_title + '.png', dpi=300)

def calculate_max_ps(ds):
    alpha = calculate_alpha(ds)     #red noise lag-one autocorrelation
    periods, pave = cal_fft(ds)
    rspec = red_noise_spectrum(len(periods), alpha)
    rspec = rspec/np.sum(rspec)

    pave_new = pave[2:]
    periods_new = periods[2:]
    rspec_new = rspec[2:]
    
    max_ind = np.argmax(pave_new)
    max_period = periods_new[max_ind]
    r = rspec_new[max_ind]
    p = pave_new[max_ind]
    sig = stats.f.cdf(p/r, dfn=2, dfd=500)
    return max_period, sig



def main():
    p_polynya = '../../SO_data/data_polynya_mean/'
    p_ice = '../../SO_data/data_siconc_w_area/'
    p_mlotst = '../../SO_data/data_mlotst/'
    p_mld = '../../SO_data/data_mld/'

    mycols = ['name', 'P_period', 'P_sig', 'C_period', 'C_sig']
    
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

        damld, dsmld = open_mld(p_mlotst, p_mld, name)
        if not check_timerange(dspolynya, damld):
            print('time range not match')
            continue

        polynya_area_ann = dsice.areacello.where(dspolynya>0).sum((dsice.areacello.dims[0],dsice.areacello.dims[1]))/1e12
        convection_area_ann = dsmld.areacello.where(damld>=2000).sum((dsmld.areacello.dims[0], dsmld.areacello.dims[1]))/1e12
        if (polynya_area_ann == 0).all():
            print("no polynya", end = '...')
            p_result = [np.nan, np.nan]
        else:
            p_results = calculate_max_ps(polynya_area_ann)
        for p in p_results:
            ps.append(p)
        if (convection_area_ann == 0).all():
            print("no convection", end = '...')
            c_result = [np.nan, np.nan]
        else:
            c_results = calculate_max_ps(convection_area_ann)
        for c in c_results:
            ps.append(c)
        # plot_spectrum(polynya_area_ann, name + '_polynya')
        # plot_spectrum(convection_area_ann, name + '_convection')
        
        df.iloc[i, 1:] = ps
        print('')
    df.to_csv('periodicity.csv', index=False)

if __name__ == "__main__":
    main()
