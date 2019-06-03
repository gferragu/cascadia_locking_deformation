#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:36:51 2019

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
from read_model import readfile
from plotting import plot_deformation_map
from plotting import plot_patches


#not sure if a meshgrid or just a 1d list of coordinates is better
lon = np.arange(-127, -118, 0.2)
lat = np.arange(40, 50, 0.2)
lonmesh, latmesh = np.meshgrid(lon,lat)
ux = np.zeros((len(lat), len(lon)))
uy = np.zeros((len(lat), len(lon)))
uz = np.zeros((len(lat), len(lon)))

myProj = Proj("+proj=utm +zone=10T, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
UTMx, UTMy = myProj(lonmesh, latmesh)
g = make_gridobj(lonmesh,latmesh, UTMx, UTMy, ux ,uy, uz)

patches_wang = readfile('wang.rupt')

patches_gaus = readfile('gaus.rupt')

#loop over every patch
#for patch p, loop over every observation point and add dispacement to the sum

from pyproj import Geod
wgs84_geod = Geod(ellps='WGS84')
for p in range(len(patches_wang.lon)):

    for j in range(len(g.utmx)):
        for i in range(len(g.utmy[0])):
            az12,az21,distx = wgs84_geod.inv(g.lon[j][i],g.lat[j][i],patches_wang.lon[p],g.lat[j][i])
            az12,az21,disty = wgs84_geod.inv(g.lon[j][i],g.lat[j][i],g.lon[j][i],patches_wang.lat[p])
            depth = patches_wang.z[p] + np.sqrt((patches_wang.ss_len[p]/2.)**2. + (patches_wang.ds_len[p]/2)**2.)
            ux,uy,uz = calc_deformation(alpha = 2./3, strike = patches_wang.strike[p], depth = depth, dip = patches_wang.dip[p],  strike_width = [-1*patches_wang.ss_len[p]/2., patches_wang.ss_len[p]/2.], dip_width = [-1*patches_wang.ds_len[p]/2., patches_wang.ds_len[p]/2.], dislocation =[patches_wang.ss_slip[p], patches_wang.ds_slip[p], 0], x = distx,y = disty)
            #next add us to g.u 
            g.ux[j][i] += ux
            g.uy[j][i] += uy
            g.uz[j][i] += uz

            
#call plotting function for map  
plot_deformation_map(lon= g.lon, lat = g.lat, uxg =  g.ux, uyg = g.uy, uzg = g.uz,title = 'gaussian model')
plt.savefig('figures/wang_fine.png')


lon = np.arange(-127, -118, 0.2)
lat = np.arange(40, 50, 0.2)
lonmesh, latmesh = np.meshgrid(lon,lat)
ux = np.zeros((len(lat), len(lon)))
uy = np.zeros((len(lat), len(lon)))
uz = np.zeros((len(lat), len(lon)))

myProj = Proj("+proj=utm +zone=10T, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
UTMx, UTMy = myProj(lonmesh, latmesh)
g = make_gridobj(lonmesh,latmesh, UTMx, UTMy, ux ,uy, uz)

for p in range(len(patches_gaus.lon)):
    for j in range(len(g.utmx)):
        for i in range(len(g.utmy[0])):
            az12,az21,distx = wgs84_geod.inv(g.lon[j][i],g.lat[j][i],patches_gaus.lon[p],g.lat[j][i])
            az12,az21,disty = wgs84_geod.inv(g.lon[j][i],g.lat[j][i],g.lon[j][i],patches_gaus.lat[p])
            depth = patches_gaus.z[p] + np.sqrt((patches_gaus.ss_len[p]/2.)**2. + (patches_gaus.ds_len[p]/2)**2.)

            ux,uy,uz = calc_deformation(alpha = 2./3, strike = patches_gaus.strike[p], depth = depth, dip = patches_gaus.dip[p],  strike_width = [-1*patches_gaus.ss_len[p]/2., patches_gaus.ss_len[p]/2.], dip_width = [-1*patches_gaus.ds_len[p]/2., patches_gaus.ds_len[p]/2.], dislocation =[patches_gaus.ss_slip[p], patches_gaus.ds_slip[p], 0],x=distx,y=disty) #x = g.utmx[j][i],y = g.utmy[j][i])
            #next add us to g.u 
#            print ux
            g.ux[j][i] += ux
            g.uy[j][i] += uy
            g.uz[j][i] += uz

print g.ux

plot_deformation_map(lon= g.lon, lat = g.lat, uxg =  g.ux, uyg = g.uy, uzg = g.uz,title = 'gaussian model')
plt.savefig('figures/gaus_fine.png')

#plot_deformation_map(longrid = g.lon.flatten(), latgrid = g.lat.flatten(), uhorizontal = 0, uvertical = 0, deformation = g.uz.flatten())

plot_patches(patches_gaus.lat, patches_gaus.lon, patches_gaus.ss_len, patches_gaus.ds_len,patches_gaus.ss_slip, label = 'strikeslip')
plt.savefig('figures/strikeslip.png')

plot_patches(patches_gaus.lat, patches_gaus.lon, patches_gaus.ss_len, patches_gaus.ds_len,patches_gaus.ds_slip, label = 'dipslip')
plt.savefig('figures/dipslip.png')







