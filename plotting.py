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
from shapely.geometry.polygon import Polygon
from matplotlib.collections import PatchCollection
from descartes import PolygonPatch

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

def plot_deformation_map(lon, lat, uxg,uyg,uzg,title, label):#and maybe arrow points of decimate later and plot
    #map of all of cascadia
    longrid = lon.flatten()
    latgrid = lat.flatten()
    ux = uxg.flatten()
    uy = uyg.flatten()
    uz = uzg.flatten()

    stamen_terrain = cimgt.StamenTerrain()#desired_tile_form="RGB")
#    stamen_terrain = cimgt.GoogleTiles(desired_tile_form='RGB', style='satellite')
    fig = plt.figure(figsize = (8,8))#6,10
    ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent([-130.01, -117.99, 40.0, 50.1])
    ax.set_extent([-128.01, -117.99, 39.8, 50.1])
    
    cmap = mpl.cm.get_cmap('viridis_r')
    normalize = mpl.colors.Normalize(vmin=-2, vmax=7)
#    normalize = mpl.colors.Normalize(vmin=-0.01, vmax=0.01)
#
    colors = [cmap(normalize(value)) for value in uz]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])

    plt.suptitle(title, fontsize = 20)

    ax.add_image(stamen_terrain, 7, alpha = 0.5)#, cmap = 'gray')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    
    xticks = [-130, -128, -126, -124, -122, -120, -118]
    yticks = [40, 42, 44, 46, 48, 50]
    ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)
    plt.scatter(longrid,latgrid, s = 30, c = colors, edgecolors = colors, transform=ccrs.PlateCarree(),alpha = 0.8,label = 'obs points')
    
    X, Y = np.meshgrid(lon[0][::2],lat[:,0][::2])
    U = uxg[::2, ::2].flatten()
    V = uyg[::2, ::2].flatten()
    q = plt.quiver(X, Y, U, V, zorder = 1000, transform=ccrs.PlateCarree())
    plt.quiverkey(q, X=0.3, Y=1.1, U=10, label='Quiver key, length = 10 m', labelpos='E')
#    X, Y = np.meshgrid(lon[0],lat[:,0])
#    U = uxg.flatten()
#    V = uxg.flatten()
#    q = plt.quiver(X, Y, U, V, zorder = 1000,transform=ccrs.PlateCarree())
#    plt.quiverkey(q, X=0.3, Y=1.1, U=5.0, label='Quiver key, length = 5 m', labelpos='E')


    fig.subplots_adjust(right=0.7)
    cbar_ax = fig.add_axes([0.75, 0.22, 0.02, 0.55])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(label, fontsize = 20)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 14)


    plt.show()

def plot_patches(lat,lon, sslen, dslen,c, label, title):#and maybe arrow points of decimate later and plot
    #map of all of cascadia
    utmxc, utmyc = myProj(lon,lat)

    lon1, lat1 = myProj(utmxc+sslen/2,utmyc+dslen/2 , inverse=True)
    lon0, lat0 = myProj(utmxc-sslen/2,utmyc-dslen/2, inverse=True)
    
    vertex_list = [[(lon0[i],lat0[i]), (lon0[i],lat1[i]), (lon1[i],lat1[i]), (lon1[i],lat0[i]), (lon0[i],lat0[i])] for i in range(len(lon1))]
    
    stamen_terrain = cimgt.StamenTerrain()#desired_tile_form="RGB")
    fig = plt.figure(figsize = (8,10))#6,10
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-128.01, -117.99, 39.8, 50.1])
    
    plt.title(title)
    
    cmap = mpl.cm.get_cmap('viridis')
    normalize = mpl.colors.Normalize(vmin=0, vmax=15)

    colors = [cmap(normalize(value)) for value in c]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    ax.add_image(stamen_terrain, 7, alpha = 0.5)#, cmap = 'gray')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    
    patches = []
    for i in range(len(vertex_list)):
        poly = Polygon(vertex_list[i])
        patches.append(PolygonPatch(poly))
    
    xticks = [-130, -128, -126, -124, -122, -120, -118]
    yticks = [40, 42, 44, 46, 48, 50]
    ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)
    
    p = PatchCollection(patches, alpha=0.4, facecolor = colors, edgecolor = colors)
    ax.add_collection(p)

    fig.subplots_adjust(right=0.7)
    cbar_ax = fig.add_axes([0.75, 0.22, 0.02, 0.55])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(label, fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 14)
    plt.show()    
    
def xsections(lon, uxwangx, uxgausx,uywangx, uygausx,uzwangx, uzgausx, title, coast):
    fig, axes = plt.subplots(3, 1, figsize = (10,10))
    axes[0].plot(lon,uxwangx, label = 'wang model')
    axes[0].plot(lon,uxgausx, label = 'gaus model')
    axes[0].set_xlabel('longitude', fontsize = 14)
    axes[0].set_ylabel('EW deformation (m)', fontsize = 14)
    axes[0].axvline(x=coast, c = 'black', label = 'coastline')
    axes[0].legend(loc=1)

    
    axes[1].plot(lon,uywangx, label = 'wang model')
    axes[1].plot(lon,uygausx, label = 'gaus model')
    axes[1].set_xlabel('longitude',fontsize = 14)
    axes[1].set_ylabel('NS deformation (m)',fontsize = 14)
    axes[1].axvline(x=coast, c = 'black', label = 'coastline')

    
    axes[2].plot(lon,uzwangx, label = 'wang model')
    axes[2].plot(lon,uzgausx, label = 'gaus model')
    axes[2].set_xlabel('longitude',fontsize = 14)
    axes[2].set_ylabel('Up down deformation (m)',fontsize = 14)
    axes[2].axvline(x=coast, c = 'black', label = 'coastline')

    plt.tight_layout()
    plt.suptitle('latitude = ' + title, fontsize = 20)
    plt.subplots_adjust(top=0.9)
    plt.show()
    plt.savefig('figures/xsection_latitude_' + title)
    
    fig, axes = plt.subplots(3, 1, figsize = (10,10))
    axes[0].plot(lon,uxwangx-uxgausx, c = 'green')
    axes[0].set_xlabel('longitude', fontsize = 14)
    axes[0].set_ylabel('EW deformation (m)', fontsize = 14)
    axes[0].axvline(x=coast, c = 'black', label = 'coastline')
    axes[0].legend(loc=1)
    
    axes[1].plot(lon,uywangx-uygausx, c = 'green')
    axes[1].set_xlabel('longitude',fontsize = 14)
    axes[1].set_ylabel('NS deformation (m)',fontsize = 14)
    axes[1].axvline(x=coast, c = 'black', label = 'coastline')
    
    axes[2].plot(lon,uzwangx-uzgausx, c = 'green')
    axes[2].set_xlabel('longitude',fontsize = 14)
    axes[2].set_ylabel('Up down deformation (m)',fontsize = 14)
    axes[2].axvline(x=coast, c = 'black', label = 'coastline')

    plt.tight_layout()
    plt.suptitle('Difference between wang and gaus latitude = ' + title, fontsize = 20)
    plt.subplots_adjust(top=0.9)
    plt.show()
    plt.savefig('figures/xsection_diff_latitude_' + title)

    
def plot_region(latlist, title):
    stamen_terrain = cimgt.StamenTerrain()#desired_tile_form="RGB")
    fig = plt.figure(figsize = (6,10))#6,10
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-128.01, -117.99, 39.8, 50.1])
    
    ax.add_image(stamen_terrain, 7, alpha = 0.5)#, cmap = 'gray')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    
    xticks = [-130, -128, -126, -124, -122, -120, -118]
    yticks = [40, 42, 44, 46, 48, 50]
    ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)
    
    lonlist = [-127. , -126.8, -126.6, -126.4, -126.2, -126. , -125.8, -125.6,
       -125.4, -125.2, -125. , -124.8, -124.6, -124.4, -124.2, -124. ,
       -123.8, -123.6, -123.4, -123.2, -123. , -122.8, -122.6, -122.4,
       -122.2, -122. , -121.8, -121.6, -121.4, -121.2, -121. , -120.8,
       -120.6, -120.4, -120.2, -120. , -119.8, -119.6, -119.4, -119.2,
       -119. , -118.8, -118.6, -118.4, -118.2]
    
    for i in range(len(latlist)):
        plt.plot(lonlist,latlist[i]*np.ones(len(lonlist)),transform=ccrs.PlateCarree())

    plt.show()    

