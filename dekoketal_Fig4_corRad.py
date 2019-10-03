#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 13:57:46 2019 (different correlation plot have same origin)

@author: Remco de Kok

Correlations between temperatures and net radiation.
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import stats
from mpl_toolkits.axes_grid1 import ImageGrid

#T @ 52m from KNMI Climate explorer
dat=Dataset('../era5_t2m.nc',mode='r')
aa=dat.variables['t2m'][:]
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
dat=Dataset('../era5_ssr.nc',mode='r')
ssr=dat.variables['ssr'][:]
dat.close()
dat=Dataset('../era5_str.nc',mode='r')
bb=ssr+dat.variables['str'][:]
dat.close()

#Create arrays for plotting
latlat=np.zeros((nx,ny))
lonlon=np.zeros((nx,ny))
dtdy=np.zeros((nt,nx,ny))

for ix in range(0,nx):
    latlat[ix,:]=lat
          
for iy in range(0,ny):
    lonlon[:,iy]=lon

#calculate total per season from monthly means for correlation, 
#weighted by number of days
seasons=['DJF','MAM','JJA','SON']
monsel=np.array([[12,1,2],[3,4,5],[6,7,8],[9,10,11]])-1
ndays=[31,28,31,30,31,30,31,31,30,31,30,31]           

cor=np.zeros((4,nx,ny))
for iseas in range(0,4):    
    pval=np.zeros((nx,ny))
    for xref in range(0,nx):
        for yref in range(0,ny):
            cora=np.zeros((nyr))
            corb=np.zeros((nyr))
            for imon in range(0,3):
                monref=monsel[iseas,imon]
                cora=cora+aa[monref::nmon,yref,xref]*ndays[monref]
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
m = Basemap(llcrnrlon=60,llcrnrlat=20,
            urcrnrlon=109,urcrnrlat=50,
            resolution='l',projection='cyl')
	   
x,y=m(lonlon[:,:],latlat[:,:])
#-------------

lonticks=np.arange(5)*10.+60
latticks=np.arange(4)*10.+20

tcol=['white','white','black','black']

fig = plt.figure(figsize=(10,8)) 
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

i=0
for ax in grid: 
    plt.sca(ax)
    ax.annotate(seasons[i],xy=(0.05,0.9),xycoords='axes fraction',size=18,color=tcol[i])

    m.drawcoastlines(linewidth=1.,color='0.') 
    im=m.contourf(x,y,cor[i,:,:],cmap='PuOr',levels=corlevs)

    for spine in plt.gca().spines.values():
        spine.set_edgecolor('0.5')
    for ip in range(0,nip):
        xx=[lonkar[ipx1[ip]],lonkar[ipx2[ip]]]
        yy=[latkar[ipy1[ip]],latkar[ipy2[ip]]]
        ax.plot(xx,yy,c='white',linewidth=2.) 

    isel=np.where((latlat>=latkar[0]) & (latlat<=latkar[1]) & (lonlon>=lonkar[0]) & (lonlon<=lonkar[1]))
    print('t2-netrad')
    corsel=cor[i,:,:]
    print('Kar '+seasons[i]+' '+str(np.mean(corsel[isel])))  


    ax.set_xticks(lonticks)
    ax.set_yticks(latticks)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(direction='in',colors='black',length=5,width=1.5,labelcolor='black',axis='both')
       
    i=i+1
    
cb = grid.cbar_axes[0].colorbar(im,ticks=[-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
cb.set_label_text('Correlation coefficient',size=15)

plt.savefig('dekoketal_corr_t2-netrad.pdf',bbox_inches='tight') 