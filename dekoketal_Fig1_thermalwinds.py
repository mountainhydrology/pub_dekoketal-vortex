#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:27:26 2019

@author: Remco de Kok

Illustrates thermal wind equation. Magnitudes are chosen for good visibility.

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

#Make a lat-lon array, and corresponding temperatures
nx=21
ny=21

lat=np.arange(ny)+20.   #Some northern mid-latitudes
lon=np.arange(nx)+10

#Coriolis parameter             
ff=2.*np.sin(lat*np.pi/180.)            

temp=np.zeros((nx,ny))
dtdy=np.zeros((nx,ny))
dtdx=np.zeros((nx,ny))
latlat=np.zeros((nx,ny))
lonlon=np.zeros((nx,ny))
for ix in range(0,nx):
    for iy in range(0,ny):
        temp[ix,iy]=280.-lat[iy]*0.5    #Temperature linearly decreasing w. lat.
        latlat[ix,iy]=lat[iy]
        lonlon[ix,iy]=lon[ix]
        
#Thermal gradients
for iy in range(0,ny):            
    dtdx[:,iy]=np.gradient(temp[:,iy])
for ix in range(0,nx):            
    dtdy[ix,:]=-np.gradient(temp[ix,:])  
    #Include Coriolis parameter from thermal wind. Other variables are constant.
    dtdy[ix,:]=dtdy[ix,:]/ff    
    dtdx[ix,:]=dtdx[ix,:]/ff

#Create good-looking temperature change
fwhm=5.
lon0=20
lat0=30
dt=np.zeros((nx,ny))
for ix in range(0,nx):
    for iy in range(0,ny):
        dt[ix,iy]=np.exp(-((lon[ix]-lon0)**4 + (lat[iy]-lat0)**4) / fwhm**4)

#Add temperature change
temp2=temp+dt*0.5   

#Calculate same thermal wind for different temperature field
dtdy2=np.zeros((nx,ny))
dtdx2=np.zeros((nx,ny))
for iy in range(0,ny):            
    dtdx2[:,iy]=np.gradient(temp2[:,iy])
for ix in range(0,nx):            
    dtdy2[ix,:]=-np.gradient(temp2[ix,:])  
    dtdy2[ix,:]=dtdy2[ix,:]/ff
    dtdx2[ix,:]=dtdx2[ix,:]/ff
         
         
#Plot State 1, State 2, and difference
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.family': 'sans-serif'})
plt.rcParams.update({'font.weight': 'light'})

fig = plt.figure(figsize=(15,5)) 
grid = ImageGrid(fig, 111,          
                 nrows_ncols=(1,3),
                 axes_pad=0.15,
                 share_all=True,
                 label_mode = "L",
                 aspect = True
                 )  


i=0
for ax in grid:
    
    plt.sca(ax)
    if i==0:
        levels=np.sort(280.-np.arange(nx))
        ax.contourf(lonlon,latlat,temp,cmap='RdBu_r',levels=levels)
        ax.quiver(lonlon[::2,::2],latlat[::2,::2],dtdy[::2,::2],dtdx[::2,::2]) 
        ax.annotate('a',xy=(0.05,0.9),xycoords='axes fraction',size=25,color='white')
        ax.set_xticks([], [])
        ax.set_yticks([], [])
        ax.set_xlabel(r'x / Longitude $\rightarrow$',size=18)
        ax.set_ylabel(r'y / Northern Latitude $\rightarrow$',size=18)
        
    if i==1:
         levels=np.sort(280.-np.arange(nx))
         ax.contourf(lonlon,latlat,temp2,cmap='RdBu_r',levels=levels)
         ax.quiver(lonlon[::2,::2],latlat[::2,::2],dtdy2[::2,::2],dtdx2[::2,::2])
         ax.annotate('b',xy=(0.05,0.9),xycoords='axes fraction',size=25,color='white')
         ax.set_xlabel(r'x / Longitude $\rightarrow$',size=18)
         
    if i==2:
        levels=np.arange(nx)/(nx-1)-0.5
        ax.contourf(lonlon,latlat,temp2-temp,cmap='RdBu_r',levels=levels)
        ax.quiver(lonlon[::2,::2],latlat[::2,::2],dtdy2[::2,::2]-dtdy[::2,::2],dtdx2[::2,::2]-dtdx[::2,::2],scale=2.)    
        ax.annotate('c',xy=(0.05,0.9),xycoords='axes fraction',size=25,color='black')
        ax.set_xlabel(r'x / Longitude $\rightarrow$',size=18)

    i=i+1    
    
plt.savefig('dekoketal_thermalwind_v2.pdf',bbox_inches='tight')  