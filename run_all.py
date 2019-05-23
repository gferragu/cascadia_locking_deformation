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

#not sure if a meshgrid or just a 1d list of coordinates is better
lon = np.arange(-127, -118, 0.5)
lat = np.arange(40, 50, 0.5)
lonmesh, latmesh = np.meshgrid(lon,lat)
ux = np.zeros((len(lon), len(lat)))
uy = np.zeros((len(lon), len(lat)))
uz = np.zeros((len(lon), len(lat)))

myProj = Proj("+proj=utm +zone=10T, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
UTMx, UTMy = myProj(lonmesh, latmesh)
g = make_gridobj(lonmesh,latmesh,UTMx, UTMy, ux,uy,uz)


patches_wang = readfile('wang.rupt')

#loop over every path
#for patch p, loop over every observation point and add dispacement to the sum

#for p in range(len(patches_wang.lon)):
    #get patch atrributes
#    
#    for i in range(len(g.lon)):
#        ux,uy,uz = calc_deformation(alpha = 2./3, strike = p.strike[p], depth = p.z[p], dip = p.dip[p],  strike_width = [utmx-ss_len/2., utmx+ss_len/2.], dip_width = [utmx-ds_len/2., utmx+ds_len/2.], dislocation, x,y)
