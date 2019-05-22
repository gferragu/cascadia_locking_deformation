#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 10:06:05 2019

@author: aklimase
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy import optimize
from okada_wrapper import dc3d0wrapper, dc3dwrapper
import matplotlib as mpl
mpl.rcParams.update({'font.size': 16})


def plot_deformation(alpha, strike, depth, dip,  strike_width, dip_width, dislocation, xmin, xmax, ymin, ymax, name):
    
    
    horizontal = []
    vert = []
    x = []
    y = []
    hx = []
    hy = []
    ux = []
    uy = []
    uz = []
    ulist = []
    
    for i in range(xmin,xmax):
        for j in range(ymin, ymax):
            strike = 0
            theta=strike-90
            theta=np.deg2rad(theta)
            R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
            xrot,yrot = R.dot(np.array([i,j]))
    
            success, u, grad_u = dc3dwrapper(alpha , [xrot, yrot, 0.0], depth, dip,  strike_width, dip_width, dislocation)                                      
            hx.append(u[0])
            hy.append(u[1])
            xrot,yrot = R.dot(np.array([i,j]))
            x.append(i)
            y.append(j)
            R2=np.array([[np.cos(-theta),-np.sin(-theta)],[np.sin(-theta),np.cos(theta)]])
            urot=R2.dot(np.array([[u[0]], [u[1]]]))
            ux.append(urot[0])
            uy.append(urot[1])
            uz.append(u[2])
            horizontal.append(float(np.sqrt(urot[0]**2 + urot[1]**2)))
            vert.append(u[2])
            ulist.append(u)
            
    
    
    #returns ux, uy, ux, at every location from xmin, xmax, ymin, ymax
    return ux,uy,uz