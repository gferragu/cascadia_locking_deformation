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
    lon	= 0 
    lat = 0
    z = 0	
    strike = 0
    dip	= 0
    risedur = 0
    ss_slip = 0	
    ds_slip = 0
    ss_len = 0	
    ds_len = 0
    rupt_time = 0
    rigidity = 0


    # The class "constructor" - It's actually an initializer 
    def __init__(self, patchnum, lon, lat, z, strike,dip,risedur,ss_slip,ds_slip,ss_len,ds_len,rupt_time,rigidity):
        self.patchnum = patchnum
        self.lon = lon
        self.lat = lat
        self.z = z
        self.strike = strike
        self.dip = dip
        self.risedur = risedur
        self.ss_slip = ss_slip
        self.ds_slip = ds_slip
        self.ss_len = ss_len
        self.ds_len = ds_len
        self.rupt_time = rupt_time
        self.rigidity = rigidity



def make_pathobj(patchnum, lon, lat, z, strike,dip,risedur,ss_slip,ds_slip,ss_len,ds_len,rupt_time,rigidity):
    pobj = path_obj(patchnum, lon, lat, z, strike,dip,risedur,ss_slip,ds_slip,ss_len,ds_len,rupt_time,rigidity)
    return pobj

#read in model file

def readfile(filename):
    data = np.genfromtxt(filename, names=True,dtype=None, delimiter="\t", encoding = None)
    pobj = make_pathobj(data['No'], data['lon'], data['lat'], data['zkm'], data['strike'],  data['dip'], data['risedura'], data['ssslipm'], data['dsslipm'],data['ss_lenm'], data['ds_lenm'], data['rupt_times'], data['rigidityPa'])
    print pobj.patchnum
    
#def define_patch():
    
    
readfile('wang.rupt')











