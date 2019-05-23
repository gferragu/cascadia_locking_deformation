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





def plot_deformation_map(longrid, latgrid, uhorizontal, uvertical):#and maybe arrow points of decimate later and plot
    #map of all of cascadia
    stamen_terrain = cimgt.StamenTerrain()#desired_tile_form="RGB")
    #stamen_terrain = cimgt.GoogleTiles(desired_tile_form='RGB', style='satellite')
    fig = plt.figure(figsize = (8,10))#6,10
    ax = plt.axes(projection=stamen_terrain.crs)
    #ax.set_extent([-130.01, -117.99, 40.0, 50.1])
    ax.set_extent([-128.01, -117.99, 40.0, 50.1])
    
    
    ax.add_image(stamen_terrain, 7, alpha = 0.5)#, cmap = 'gray')
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.tick_params(axis='both', which='both', length = 5, width = 1)
    
    xticks = [-130, -128, -126, -124, -122, -120, -118]
    yticks = [40, 42, 44, 46, 48, 50]
    ax.gridlines(xlocs=xticks, ylocs=yticks, draw_labels = True)
#    plt.scatter(g.lon,g.lat, s= 20, c = 'black', transform=ccrs.Geodetic(),label = 'obs points')
    

plot_deformation_map()
    
    
#def compare_models():