#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:36:50 2019

@author: aklimase
"""
from okada import calc_deformation
import cartopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from read_model import make_gridobj
from pyproj import Proj

#not sure if a meshgrid or just a 1d list of coordinates is better
lon = np.arange(-127, -120, 0.5)
lat = np.arange(40, 50, 0.5)
lonmesh, latmesh = np.meshgrid(lon,lat)
ux = np.zeros((len(lat), len(lon)))
uy = np.zeros((len(lat), len(lon)))
uz = np.zeros((len(lat), len(lon)))

myProj = Proj("+proj=utm +zone=10T, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
UTMx, UTMy = myProj(lonmesh, latmesh)
g = make_gridobj(lonmesh,latmesh,UTMx, UTMy, ux,uy,uz)

def plot_deformation_map(lon, lat, ux,uy,uz):#and maybe arrow points of decimate later and plot
    #map of all of cascadia
    longrid = lon.flatten()
    latgrid = lat.flatten()
    
    
    stamen_terrain = cimgt.StamenTerrain()#desired_tile_form="RGB")
#    stamen_terrain = cimgt.GoogleTiles(desired_tile_form='RGB', style='satellite')
    fig = plt.figure(figsize = (8,10))#6,10
    ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent([-130.01, -117.99, 40.0, 50.1])
    ax.set_extent([-128.01, -117.99, 40.0, 50.1])
    
    cmap = mpl.cm.get_cmap('viridis')
    normalize = mpl.colors.Normalize(vmin=min(uz), vmax=max(uz))
#    normalize = mpl.colors.Normalize(vmin=-0.01, vmax=0.01)
#
    colors = [cmap(normalize(value)) for value in uz]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])

    
    
    ax.add_image(stamen_terrain, 7, alpha = 0.5)#, cmap = 'gray')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    
#
#    horiz=(ux**2+uy**2)**0.5
#    plt.quiver(longrid,latgrid,ux,uy,color='k')
#
    
    xticks = [-130, -128, -126, -124, -122, -120, -118]
    yticks = [40, 42, 44, 46, 48, 50]
    ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)
    plt.scatter(longrid,latgrid, s= 30, c = colors, edgecolors = None, transform=ccrs.PlateCarree(),alpha = 0.7,label = 'obs points')
    
    X, Y = np.meshgrid(lon[0],lat[:,0])
    U = ux
    V = uy
    q = plt.quiver(X,Y, U, V, zorder = 1000,transform=ccrs.PlateCarree())
    plt.quiverkey(q, X=0.3, Y=1.1, U=1, label='Quiver key', labelpos='E')

    
    
    fig.subplots_adjust(right=0.7)
    cbar_ax = fig.add_axes([0.75, 0.12, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"deformation (m)", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 14)



    plt.show()

#plot_deformation_map(g.lon, g.lat, [0],[0])
    
    
#def compare_models():


