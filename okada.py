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


def calc_deformation(alpha, strike, depth, dip,  strike_width, dip_width, dislocation, x,y):
    

    theta=strike-90
    theta=np.deg2rad(theta)
    R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    R2=np.array([[np.cos(-theta),-np.sin(-theta)],[np.sin(-theta),np.cos(theta)]])

    xrot,yrot = R.dot(x,y)

    success, u, grad_u = dc3dwrapper(alpha , [xrot, yrot, 0.0], depth, dip,  strike_width, dip_width, dislocation)                                      

    urot=R2.dot(np.array([[u[0]], [u[1]]]))
    u[0]=urot[0]
    u[1]=urot[1]
    
    ux = u[0]
    uy = u[1]
    uz = u[2]
    
    
    #returns ux, uy, ux, at every location from xmin, xmax, ymin, ymax
    return ux,uy,uz