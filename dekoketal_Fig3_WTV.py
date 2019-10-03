#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 13:57:46 2019 (different correlation plot have same origin)

@author: Remco de Kok

Correlations between mean temperatures and wind speeds, and temperatures and zonal shear.
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import stats
from mpl_toolkits.axes_grid1 import ImageGrid

#U @ 300 hPa from KNMI Climate explorer
dat=Dataset('../era5_u300.nc',mode='r')
aa=dat.variables['u'][:]
lon=dat.variables['lon'][:]
lat=dat.variables['lat'][:]
dat.close()
temp=aa
aa=aa[:,0,:,:]

nx=len(lon)
ny=len(lat)
dims=np.shape(aa)
nt=dims[0]
nmon=12
nyr=int(nt/nmon)

#T @ 2m from KNMI Climate explorer
dat=Dataset('../era5_t2m.nc',mode='r')
bb=dat.variables['t2m'][:]
dat.close()

#Arrays for plotting and selecting
latlat=np.zeros((nx,ny))
lonlon=np.zeros((nx,ny))

for ix in range(0,nx):
    latlat[ix,:]=lat
          
for iy in range(0,ny):
    lonlon[:,iy]=lon
          
latlat2=np.transpose(latlat) 
lonlon2=np.transpose(lonlon)  

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
                   
isel1=np.where((latlat2>=latbox1[0]) & (latlat2<=latbox1[1]) & (lonlon2>=lonbox1[0]) & (lonlon2<=lonbox1[1]))
isel2=np.where((latlat2>=latbox2[0]) & (latlat2<=latbox2[1]) & (lonlon2>=lonbox2[0]) & (lonlon2<=lonbox2[1]))
iselkar=np.where((latlat2>=latkar[0]) & (latlat2<=latkar[1]) & (lonlon2>=lonkar[0]) & (lonlon2<=lonkar[1]))

#Shifted reference boxes
shiftlat=15.
shiftlon=36.
isel1a=np.where((latlat2>=latbox1[0]+shiftlat) & (latlat2<=latbox1[1]+shiftlat) & (lonlon2>=lonbox1[0]+shiftlon) & (lonlon2<=lonbox1[1]+shiftlon))
isel2a=np.where((latlat2>=latbox2[0]+shiftlat) & (latlat2<=latbox2[1]+shiftlat) & (lonlon2>=lonbox2[0]+shiftlon) & (lonlon2<=lonbox2[1]+shiftlon))
iselkara=np.where((latlat2>=latkar[0]+shiftlat) & (latlat2<=latkar[1]+shiftlat) & (lonlon2>=lonkar[0]+shiftlon) & (lonlon2<=lonkar[1]+shiftlon))
    
#Calculate for seasons. Here, we only took JJA
seasons=['DJF','MAM','JJA','SON']
monsel=np.array([[12,1,2],[3,4,5],[6,7,8],[9,10,11]])-1
ndays=[31,28,31,30,31,30,31,31,30,31,30,31]               

for iseas in range(2,3):    
    cor_figa=np.zeros((nx,ny))   #correlation coefficients
    cor_figb=np.zeros((nx,ny))  
    cor_figc=np.zeros((nx,ny))   
    cor_figd=np.zeros((nx,ny))   
    
    cora=np.zeros((nyr,ny,nx))  #arrays for correlation
    corb=np.zeros((nyr,ny,nx))
    kzi=np.zeros((nyr))     #zonal shear
    tkar=np.zeros((nyr))    #mean in reference box
    kzishift=np.zeros((nyr))
    tkarshift=np.zeros((nyr))    
    for imon in range(0,3):
        monref=monsel[iseas,imon]
        cora=cora+aa[monref::nmon,:,:]*ndays[monref]
        corb=corb+bb[monref::nmon,:,:]*ndays[monref]
        for iyr in range(0,nyr):
            #Calculate indices and means per year
            corasel=cora[iyr,:,:] 
            corbsel=corb[iyr,:,:]   
            kzi[iyr]=np.mean(corasel[isel1])-np.mean(corasel[isel2])
            tkar[iyr]=np.mean(corbsel[iselkar])
            kzishift[iyr]=np.mean(corasel[isel1a])-np.mean(corasel[isel2a])
            tkarshift[iyr]=np.mean(corbsel[iselkara])
            
    for xref in range(0,nx):
        for yref in range(0,ny):
            cora2=np.zeros((nyr))
            corb2=np.zeros((nyr))
            
            for imon in range(0,3):
                monref=monsel[iseas,imon]
                cora2=cora2+aa[monref::nmon,yref,xref]*ndays[monref]
                corb2=corb2+bb[monref::nmon,yref,xref]*ndays[monref]
                                
            r,p=stats.spearmanr(cora2,tkar)
            cor_figa[xref,yref]=r

            r,p=stats.spearmanr(cora2,tkarshift)
            cor_figc[xref,yref]=r

            r,p=stats.spearmanr(corb2,kzi)
            cor_figb[xref,yref]=r                    
    
            r,p=stats.spearmanr(corb2,kzishift)
            cor_figd[xref,yref]=r   
  
#Plot everything and print mean values
                  
m = Basemap(llcrnrlon=30,llcrnrlat=10,
            urcrnrlon=140,urcrnrlat=70,
            resolution='l',projection='cyl')

	   
x,y=m(lonlon[:,:],latlat[:,:])

nlevs=10
corlevs=np.arange(nlevs+1)*2/(nlevs)-1
                 
lonticks=np.arange(8)*15.+30
latticks=np.arange(4)*15.+15
#-------------

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

i=0
for ax in grid:
    plt.sca(ax)


    
    if i==0:
        m.drawcoastlines(linewidth=1.,color='0.') 
        im=m.contourf(x,y,cor_figa,cmap='PuOr',levels=corlevs)
        ax.annotate('a',xy=(0.93,0.05),xycoords='axes fraction',size=22,color='black')

        for spine in plt.gca().spines.values():
            spine.set_edgecolor('0.5')
        for ip in range(0,nip):
            xx=[lonbox1[ipx1[ip]],lonbox1[ipx2[ip]]]
            yy=[latbox1[ipy1[ip]],latbox1[ipy2[ip]]]
            ax.plot(xx,yy,c='white',linewidth=2.)  
            xx=[lonbox2[ipx1[ip]],lonbox2[ipx2[ip]]]
            yy=[latbox2[ipy1[ip]],latbox2[ipy2[ip]]]
            ax.plot(xx,yy,'white',linewidth=2.)  

        print('t2-U300')
        isel=np.where((latlat>=latbox1[0]) & (latlat<=latbox1[1]) & (lonlon>=lonbox1[0]) & (lonlon<=lonbox1[1]))
        print('BOX 1 '+seasons[i]+' '+str(np.mean(cor_figa[isel])))
        isel=np.where((latlat>=latbox2[0]) & (latlat<=latbox2[1]) & (lonlon>=lonbox2[0]) & (lonlon<=lonbox2[1]))
        print('BOX 2 '+seasons[i]+' '+str(np.mean(cor_figa[isel])))    
    
    if i==1:
        m.drawcoastlines(linewidth=1.,color='0.') 
        im=m.contourf(x,y,cor_figb,cmap='PuOr',levels=corlevs)
        ax.annotate('b',xy=(0.93,0.05),xycoords='axes fraction',size=22,color='black')
        
        for spine in plt.gca().spines.values():
            spine.set_edgecolor('0.5')
        for ip in range(0,nip):
            xx=[lonkar[ipx1[ip]],lonkar[ipx2[ip]]]
            yy=[latkar[ipy1[ip]],latkar[ipy2[ip]]]
            ax.plot(xx,yy,c='white',linewidth=2.) 
        isel=np.where((latlat>=latkar[0]) & (latlat<=latkar[1]) & (lonlon>=lonkar[0]) & (lonlon<=lonkar[1]))
        print('t2-DelU')
        print('Kar '+seasons[i]+' '+str(np.mean(cor_figb[isel])))              
        
    if i==2:
        m.drawcoastlines(linewidth=1.,color='0.') 
        im=m.contourf(x,y,cor_figc,cmap='PuOr',levels=corlevs)
        ax.annotate('c',xy=(0.93,0.05),xycoords='axes fraction',size=22,color='black')        
        for spine in plt.gca().spines.values():
            spine.set_edgecolor('0.5')
        for ip in range(0,nip):
            xx=np.array([lonbox1[ipx1[ip]],lonbox1[ipx2[ip]]])+36.
            yy=np.array([latbox1[ipy1[ip]],latbox1[ipy2[ip]]])+15.
            ax.plot(xx,yy,c='white',linewidth=2.)  
            xx=np.array([lonbox2[ipx1[ip]],lonbox2[ipx2[ip]]])+36.
            yy=np.array([latbox2[ipy1[ip]],latbox2[ipy2[ip]]])+15.
            ax.plot(xx,yy,'white',linewidth=2.)  
        print('t2m-U300')
        isel=np.where((latlat>=latbox1[0]+15.) & (latlat<=latbox1[1]+15) & (lonlon>=lonbox1[0]+36) & (lonlon<=lonbox1[1]+36))
        print('BOX 1 shifted'+seasons[i]+' '+str(np.mean(cor_figc[isel])))
        isel=np.where((latlat>=latbox2[0]+15.) & (latlat<=latbox2[1]+15) & (lonlon>=lonbox2[0]+36) & (lonlon<=lonbox2[1]+36))
        print('BOX 2 shifted'+seasons[i]+' '+str(np.mean(cor_figc[isel])))   

    
    if i==3:
        m.drawcoastlines(linewidth=1.,color='0.') 
        im=m.contourf(x,y,cor_figd,cmap='PuOr',levels=corlevs)
        ax.annotate('d',xy=(0.93,0.05),xycoords='axes fraction',size=22,color='black')
        for spine in plt.gca().spines.values():
            spine.set_edgecolor('0.5')
        for ip in range(0,nip):
            xx=np.array([lonkar[ipx1[ip]],lonkar[ipx2[ip]]])+36.
            yy=np.array([latkar[ipy1[ip]],latkar[ipy2[ip]]])+15.
            ax.plot(xx,yy,c='white',linewidth=2.) 
        isel=np.where((latlat>=latkar[0]+15) & (latlat<=latkar[1]+15) & (lonlon>=lonkar[0]+36) & (lonlon<=lonkar[1]+36))
        print('t2-DelU')
        print('Kar shifted'+seasons[i]+' '+str(np.mean(cor_figd[isel])))              

    ax.set_xticks(lonticks)
    ax.set_yticks(latticks)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(direction='in',colors='black',length=5,width=1.5,labelcolor='black',axis='both')
        
            
    i=i+1

cb = grid.cbar_axes[0].colorbar(im,ticks=[-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1])
cb.set_label_text('Correlation coefficient',size=15)

plt.savefig('dekoketal_vortices.pdf',bbox_inches='tight')                   