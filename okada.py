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


alpha = 0.6
depth = 0.0
dip = 65
strike_width = [-0.0,100]
dip_width = [-15, 0.0]
dislocation = [5., 0.0, 0.0]
xmin = -50
xmax = 50
ymin = -25
ymax = 125
name = 'NSstrike'

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
            
    horizontalq = []
    vertq = []
    xq = []
    yq = []
    hxq = []
    hyq = []
    for i in range(xmin,xmax,10):
        for j in range(ymin, ymax, 10):
            strike = 0
            theta=strike-90
            theta=np.deg2rad(theta)
            R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
            xrot,yrot = R.dot(np.array([i,j]))
    
            success, u, grad_u = dc3dwrapper(alpha , [xrot, yrot, 0.0], depth, dip,  strike_width, dip_width, dislocation)                                      
                                    
            hxq.append(u[0])
            hyq.append(u[1])
            xq.append(i)
            yq.append(j)
            R2=np.array([[np.cos(-theta),-np.sin(-theta)],[np.sin(-theta),np.cos(theta)]])
            urot=R2.dot(np.array([[u[0]], [u[1]]]))
            horizontalq.append(float(np.sqrt(urot[0]**2 + urot[1]**2)))
            vertq.append(u[2])
    
    
    
    #horizontal = np.sqrt(np.pow(ux,2.) + np.pow(uy,2.))      
    cmap = mpl.cm.get_cmap('plasma')
    normalize = mpl.colors.Normalize(vmin=min(horizontal),vmax = max(horizontal))
    colors = [cmap(normalize(value)) for value in horizontal]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    
    fig= plt.figure(figsize = (10,8))
    plt.title(r'horizontal displacement')
    plt.ylabel(r'y (km)')
    plt.xlabel(r'x (km)')
    
    plt.scatter(x,y, s = 12,c = colors, edgecolors = colors, marker = 'o', alpha = 0.6)
    fig.subplots_adjust(right=0.85)
    
    X, Y = np.meshgrid(range(xmin, xmax,5),range(ymin, ymax,5))
    U = hxq
    V = hyq
    q = plt.quiver(xq, yq, U, V)
    plt.quiverkey(q, X=0.3, Y=1.1, U=1, label='Quiver key, length = 1', labelpos='E')
    
    
    
    cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"displacement (m)", fontsize = 22)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    plt.show()
    plt.savefig(workingdir + 'horizontal_NS.png')
    
    #vert = uy
    cmap = mpl.cm.get_cmap('plasma')
    normalize = mpl.colors.Normalize(vmin=min(vert),vmax = max(vert))
    colors = [cmap(normalize(value)) for value in vert]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    
    
    fig = plt.figure(figsize = (10,8))
    plt.title(r'vertical displacement')
    plt.ylabel(r'y (km)')
    plt.xlabel(r'x (km)')
    plt.scatter(x,y, s = 12,c = colors, edgecolors = colors, alpha = 0.6)
    fig.subplots_adjust(right=0.85)
    
    X, Y = np.meshgrid(range(xmin, xmax,5),range(ymin, ymax,5))
    U = hxq
    V = hyq
    q = plt.quiver(xq, yq, U, V)
    plt.quiverkey(q, X=0.3, Y=1.1, U=1, label='Quiver key, length = 1', labelpos='E')
    
    
    cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"displacement(m)", fontsize = 18)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    plt.show()
    plt.savefig(workingdir + 'vert_NS.png')
    
    
    
    
    
    #now dot with lOS
    look = [-0.53, -0.11, 0.84]
    unwrapped = []
    for i in range(len(ulist)):
        unwrapped.append(np.dot(ulist[i],look))
    #    print unwrapped
    
    
    cmap = mpl.cm.get_cmap('viridis')
    normalize = mpl.colors.Normalize(vmin=min(unwrapped),vmax = max(unwrapped))
    colors = [cmap(normalize(value)) for value in unwrapped]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    
    fig = plt.figure(figsize = (10,8))
    plt.title(r'LOS displacement')
    plt.ylabel(r'y (km)')
    plt.xlabel(r'x (km)')
    plt.scatter(x,y, s = 12,c = colors, edgecolors = colors, alpha = 0.6)
    fig.subplots_adjust(right=0.85)
    
    cbar_ax = fig.add_axes([0.86, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"LOS displacement(m)", fontsize = 18)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    plt.show()
    plt.savefig(workingdir + 'LOS_unwrapped_NS.png')
    
    wrapped = []
    for i in range(len(unwrapped)):
        for n in range(0,30):
            if (n*0.163)<= abs(unwrapped[i])<(n+1)*0.163:
                wrapped.append((abs(unwrapped[i]) - n*0.163))
    
    
    
    cmap = mpl.cm.get_cmap('gist_rainbow')
    normalize = mpl.colors.Normalize(vmin=min(wrapped),vmax = max(wrapped))
    colors = [cmap(normalize(value)) for value in wrapped]
    s_m = mpl.cm.ScalarMappable(cmap = cmap, norm=normalize)
    s_m.set_array([])
    
    
    fig = plt.figure(figsize = (10,8))
    plt.title(r'wrapped LOS displacement')
    plt.ylabel(r'y (km)')
    plt.xlabel(r'x (km)')
    plt.scatter(x,y, s = 12,c = colors, edgecolors = colors, alpha = 0.6)
    fig.subplots_adjust(right=0.85)
    
    cbar_ax = fig.add_axes([0.86, 0.15, 0.02, 0.75])
    cbar = plt.colorbar(s_m, cax=cbar_ax)
    cbar.set_label(ur"LOS displacement(m)", fontsize = 18)#"$azimuth$ (\u00b0)"
    cbar.ax.tick_params(labelsize = 18)
    plt.show()
    plt.savefig(workingdir + 'LOS_wrapped_NS.png')
    
    
    plt.close('all')

