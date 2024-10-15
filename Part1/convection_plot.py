import numpy as np
import xarray as xr
import os
import pickle
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.path as mpath

def plottest(pf, circle, savef):
    with open(pf, 'rb') as f:
        data_dict = pickle.load(f)
    pltm = xr.where(data_dict['mldmax']>=1000, 1, np.nan)
    pltf = data_dict['mld2kfq']
    pltx = data_dict['pltx']
    plty = data_dict['plty']

    fig = plt.figure()
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180, 180, -90, -50], ccrs.PlateCarree()) 
    
    try:
        im0 = ax.pcolormesh(pltx, plty, pltm, vmin = 0, vmax=1.2, 
                            transform=ccrs.PlateCarree(), cmap=plt.cm.gray)
        ax.add_feature(cfeature.LAND, zorder=1)
        ax.add_feature(cfeature.OCEAN)
        ax.set_boundary(circle, transform=ax.transAxes)
        im = ax.pcolormesh(pltx, plty, pltf, vmin = 0, vmax=0.025, 
                           transform=ccrs.PlateCarree(), cmap=plt.cm.Spectral_r)
        cbar = plt.colorbar(im) 
        fig.savefig(savef)
        return fig
    except Exception as e:
        return e

def main():
    path = "data_mld/"
    # cwd = os.getcwd()

    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    for f in os.listdir(path):
        opf = path + '/' + f
        fname = os.path.splitext(f)
        savef = 'fig_mld_test/' + fname[0] + '.png'
        # print(savef)
        # print(opf)
        print(plottest(opf, circle, savef))
        # break

if __name__ == '__main__':
    main()

# fig = plt.figure()
# fig.savefig('test.png')