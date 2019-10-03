#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 13:57:46 2019 (different correlation plot have same origin)

@author: Remco de Kok

Correlations between thermal gradient and wind speeds.
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import stats
from mpl_toolkits.axes_grid1 import ImageGrid

#T @ 500 hPa from KNMI Climate explorer
dat=Dataset('../era5_t500.nc',mode='r')
aa=dat.variables['t'][:]
lon=dat.variables['lon'][:]
lat=dat.variables['lat'][:]
dat.close()

nx=len(lon)
ny=len(lat)
dims=np.shape(aa)
nt=dims[0]
nmon=12
nyr=int(nt/nmon)

#U @ 300 hPa from KNMI Climate explorer
dat=Dataset('../era5_u300.nc',mode='r')
bb=dat.variables['u'][:]
dat.close()
bb=bb[:,0,:,:]  #bb = variable 2, cut our height dimension in array

#Create arrays for plotting
latlat=np.zeros((nx,ny))
lonlon=np.zeros((nx,ny))
dtdy=np.zeros((nt,nx,ny))

for ix in range(0,nx):
    latlat[ix,:]=lat
    for it in range(0,nt):
        #Zonal gradient. 
        #Conversion from lat to length is constant, so does not affect correlation
        dtdy[it,ix,:]=np.gradient(aa[it,0,:,ix])
          
for iy in range(0,ny):
    lonlon[:,iy]=lon

#calculate total per season from monthly means for correlation, 
#weighted by number of days
seasons=['DJF','MAM','JJA','SON']
monsel=np.array([[12,1,2],[3,4,5],[6,7,8],[9,10,11]])-1
ndays=[31,28,31,30,31,30,31,31,30,31,30,31]       

aa=dtdy     #Variable 1        

cor=np.zeros((4,nx,ny))
for iseas in range(0,4):    
    pval=np.zeros((nx,ny))
    for xref in range(0,nx):
        for yref in range(0,ny):
            cora=np.zeros((nyr))
            corb=np.zeros((nyr))
            for imon in range(0,3):
                monref=monsel[iseas,imon]
                cora=cora+aa[monref::nmon,xref,yref]*ndays[monref]
                corb=corb+bb[monref::nmon,yref,xref]*ndays[monref]

            #Do correlation for each point
            r,p=stats.spearmanr(cora,corb)
            cor[iseas,xref,yref]=r
            pval[xref,yref]=p
    
    
nlevs=10
corlevs=np.arange(nlevs+1)*2/(nlevs)-1
#-------------

#karakoram vortex reference boxes
latbox1=[40,50]
lonbox1=[52.5,86.25]
latbox2=[20,32.5]
lonbox2=[52.5,93.75]

latkar=[35.,36.]
lonkar=[74.,76.]


ipx1=[0,0,1,1]
ipx2=[1,0,0,1]
ipy1=[0,0,1,0]
ipy2=[0,1,1,1]
nip=len(ipx1)

#Plot

plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'font.family': 'sans-serif'})
plt.rcParams.update({'font.weight': 'light'})

#Plotting range and projection
m = Basemap(llcrnrlon=0,llcrnrlat=-80,
            urcrnrlon=185,urcrnrlat=80,
            resolution='l',projection='cyl')

#Have consistent longitude definition
lonlon[lonlon<0]=lonlon[lonlon<0]+360.
	   
x,y=m(lonlon[:,:],latlat[:,:])

fig = plt.figure(figsize=(10,10)) 
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(2,2),
                 axes_pad=0.15,
                 share_all=True,
                 cbar_location="bottom",
                 cbar_mode="single",
                 cbar_size="3%",
                 cbar_pad=0.15,
                 label_mode = "L",
                 )  

seasons=['DJF','MAM','JJA','SON']
nseas=len(seasons)

lonticks=np.arange(7)*30.
latticks=np.arange(5)*30.-60

i=0
for ax in grid:
    
    plt.sca(ax)
    m.drawcoastlines(linewidth=1.,color='0.') 
    im=m.contourf(x,y,cor[i,:,:],cmap='PuOr',levels=corlevs)
    ax.annotate(seasons[i],xy=(0.05,0.9),xycoords='axes fraction',size=18,color='white')

    ax.set_xticks(lonticks)
    ax.set_yticks(latticks)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(direction='in',colors='white',length=5,width=1.5,labelcolor='black',axis='both')
    
    for spine in plt.gca().spines.values():
        spine.set_edgecolor('0.5')
    for ip in range(0,nip):
        xx=[lonbox1[ipx1[ip]],lonbox1[ipx2[ip]]]
        yy=[latbox1[ipy1[ip]],latbox1[ipy2[ip]]]
        ax.plot(xx,yy,c='white',linewidth=1.5)  
        xx=[lonbox2[ipx1[ip]],lonbox2[ipx2[ip]]]
        yy=[latbox2[ipy1[ip]],latbox2[ipy2[ip]]]
        ax.plot(xx,yy,'white',linewidth=1.5)  

    #For values in Table S1
    print('u300-dtdy500')
    isel=np.where((latlat>=latbox1[0]) & (latlat<=latbox1[1]) & (lonlon>=lonbox1[0]) & (lonlon<=lonbox1[1]))
    corsel=cor[i,:,:]
    print('BOX 1 '+seasons[i]+' '+str(np.mean(corsel[isel])))
    isel=np.where((latlat>=latbox2[0]) & (latlat<=latbox2[1]) & (lonlon>=lonbox2[0]) & (lonlon<=lonbox2[1]))
    print('BOX 2 '+seasons[i]+' '+str(np.mean(corsel[isel])))    
    
    i=i+1
    

cb = grid.cbar_axes[0].colorbar(im,ticks=[-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
cb.set_label_text('Correlation coefficient',size=15)

plt.savefig('dekoketal_corr_u300-dtdy500.pdf',bbox_inches='tight') 