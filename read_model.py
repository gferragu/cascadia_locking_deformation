#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 11:24:39 2019

@author: aklimase
"""
import numpy as np
from pyproj import Proj

class path_obj(object):
    patchnum = 0
    lon= 0 
    lat = 0
    utmx = 0
    utmy = 0
    z = 0	
    strike = 0
    dip	= 0
    rise = 0
    dur = 0
    ss_slip = 0	
    ds_slip = 0
    ss_len = 0	
    ds_len = 0
    rupt_time = 0
    rigidity = 0


    # The class "constructor" - It's actually an initializer 
    def __init__(self, patchnum, lon, lat, utmx, utmy, z, strike,dip,rise, dur,ss_slip,ds_slip,ss_len,ds_len,rupt_time,rigidity):
        self.patchnum = patchnum
        self.lon = lon
        self.lat = lat
        self.utmx = utmx
        self.utmy = utmy
        self.z = z
        self.strike = strike
        self.dip = dip
        self.rise = rise
        self.dur = dur
        self.ss_slip = ss_slip
        self.ds_slip = ds_slip
        self.ss_len = ss_len
        self.ds_len = ds_len
        self.rupt_time = rupt_time
        self.rigidity = rigidity



def make_pathobj(patchnum, lon, lat, utmx, utmy, z, strike,dip,rise,dur,ss_slip,ds_slip,ss_len,ds_len,rupt_time,rigidity):
    pobj = path_obj(patchnum, lon, lat, utmx, utmy, z, strike,dip,rise,dur,ss_slip,ds_slip,ss_len,ds_len,rupt_time,rigidity)
    return pobj

#read in model file

def readfile(filename):
    data = np.genfromtxt(filename, names=True,dtype=None, delimiter=None, encoding = None)
    myProj = Proj("+proj=utm +zone=10T, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    UTMx, UTMy = myProj(data['lon'], data['lat'])
    pobj = make_pathobj(data['No'], data['lon'], data['lat'], UTMx, UTMy, data['zkm'], data['strike'],  data['dip'], data['rise'], data['dura'], data['ssslipm'], data['dsslipm'],data['ss_lenm'], data['ds_lenm'], data['rupt_times'], data['rigidityPa'])
    return pobj
    

    
patches_wang = readfile('wang.rupt')
patches_gaus = readfile('gaus.rupt')






class outgrid_obj(object):
    lon= 0 
    lat = 0
    utmx = 0
    utmy = 0
    ux = 0
    uy = 0
    uz = 0


    # The class "constructor" - It's actually an initializer 
    def __init__(self, lon, lat, utmx, utmy, ux,uy,uz):
        self.lon = lon
        self.lat = lat
        self.utmx = utmx
        self.utmy = utmy
        self.ux = ux
        self.uy = uy
        self.uz = uz


def make_gridobj(lon, lat, utmx, utmy, ux,uy,uz):
    gobj = outgrid_obj(lon, lat, utmx, utmy, ux,uy,uz)
    return gobj

#not sure if a meshgrid or just a 1d list of coordinates is better
lon = np.arange(-127, -120, 1)
lat = np.arange(40, 50, 1)
lonmesh, latmesh = np.meshgrid(lon,lat)
ux = np.zeros((len(lat), len(lon)))
uy = np.zeros((len(lat), len(lon)))
uz = np.zeros((len(lat), len(lon)))

myProj = Proj("+proj=utm +zone=10T, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
UTMx, UTMy = myProj(lonmesh, latmesh)
g = make_gridobj(lonmesh,latmesh,UTMx, UTMy, ux,uy,uz)


