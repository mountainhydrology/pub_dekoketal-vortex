#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 13:08:52 2019

@author: Remco de Kok

Calculates thermodynamic components of heating for two years,
and looks at differences.
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

pref=1000.  #mbar
kk=0.286     #R/c_p

#Hourly data from ERA5 
dat=Dataset('../Karakoram_hourly/1987.nc',mode='r')
u1=dat.variables['u'][:]
v1=dat.variables['v'][:]
w1=dat.variables['w'][:] * 0.01 #Pa to mbar
t1=dat.variables['t'][:]
pres=dat.variables['level'][:]
lat=dat.variables['latitude'][:]
lon=dat.variables['longitude'][:]
dat.close()

dat=Dataset('../Karakoram_hourly/1994.nc',mode='r')
u2=dat.variables['u'][:]
v2=dat.variables['v'][:]
w2=dat.variables['w'][:] * 0.01 #Pa to mbar
t2=dat.variables['t'][:]
dat.close()

dims=np.shape(t1)
nt=dims[0]
npres=dims[1]
nlat=dims[2]
nlon=dims[3]

#convert angles to distance, assuming a sphere
distlat=(lat-lat[0])*111.e3
distlon=(lon-lon[0])*111.e3 * np.cos(lat*np.pi/180.)    

#Factors for the potential temperature
pfac=(pref/pres)**kk
pott1=t1*pfac
pott2=t2*pfac

#Pressures in ERA5 data
pres=np.array([200,225,250,300,350,400,450,500,550])
preslabel=['500','400','300','200']
hgt=np.zeros((npres))
for ipres in range(npres-2,-1,-1):
    hgt[ipres]=hgt[ipres+1]+np.log(pres[ipres+1]-pres[ipres])

#centre point     
ix=4
iy=4

dt0=t2[0,:,ix,iy]-t1[0,:,ix,iy]     #Initial temp. difference.

      
#Calculate all thermodynamic components   
#Mean temperatures and winds   
mean_t1=np.mean(t1[:,:,ix,iy],axis=0)
mean_t2=np.mean(t2[:,:,ix,iy],axis=0)
mean_u1=np.mean(u1[:,:,ix,iy],axis=0)
mean_u2=np.mean(u2[:,:,ix,iy],axis=0)
mean_v1=np.mean(v1[:,:,ix,iy],axis=0)
mean_v2=np.mean(v2[:,:,ix,iy],axis=0)
mean_w1=np.mean(w1[:,:,ix,iy],axis=0)
mean_w2=np.mean(w2[:,:,ix,iy],axis=0)

#Heating rates
dt1=np.zeros((nt,npres))
dpottdp1=np.zeros((nt,npres))
adh1=np.zeros((nt,npres))
htad1=np.zeros((nt,npres))
dh1a=np.zeros((nt,npres))
dh1b=np.zeros((nt,npres))

dt2=np.zeros((nt,npres))
dpottdp2=np.zeros((nt,npres))
adh2=np.zeros((nt,npres))
htad2=np.zeros((nt,npres))
dh2a=np.zeros((nt,npres))
dh2b=np.zeros((nt,npres))

#Integrated temperatures for different components
temp_adh1=np.zeros((nt,npres))
temp_htad1=np.zeros((nt,npres))
temp_dh1a=np.zeros((nt,npres))
temp_dh1b=np.zeros((nt,npres))

temp_adh2=np.zeros((nt,npres))
temp_htad2=np.zeros((nt,npres))
temp_dh2a=np.zeros((nt,npres))
temp_dh2b=np.zeros((nt,npres))

for it in range(0,nt):
    
    dpottdp1[it,:]=np.gradient(pott1[it,:,ix,iy],pres)
    adh1[it,:]=w1[it,:,ix,iy]*dpottdp1[it,:]/pfac

    dpottdp2[it,:]=np.gradient(pott2[it,:,ix,iy],pres)
    adh2[it,:]=w2[it,:,ix,iy]*dpottdp2[it,:]/pfac
        
    if it >0:
        dt1[it,:]=(t1[it,:,ix,iy]-t1[it-1,:,ix,iy])/3600. #per second
        dt2[it,:]=(t2[it,:,ix,iy]-t2[it-1,:,ix,iy])/3600. #per second
        
    for ipres in range(0,npres):
        dtdx1=np.gradient(t1[it,ipres,:,iy],distlat)
        dtdy1=np.gradient(t1[it,ipres,ix,:],distlon)
        htad1[it,ipres]=(dtdx1[ix]*v1[it,ipres,ix,iy])+(dtdy1[ix]*u1[it,ipres,ix,iy])
        
        dudx1=np.gradient(u1[it,ipres,:,iy],distlat)
        dvdy1=np.gradient(v1[it,ipres,ix,:],distlon)
             
        dtdx2=np.gradient(t2[it,ipres,:,iy],distlat)
        dtdy2=np.gradient(t2[it,ipres,ix,:],distlon)
        htad2[it,ipres]=(dtdx2[ix]*v2[it,ipres,ix,iy])+(dtdy2[ix]*u2[it,ipres,ix,iy])
        
        dudx2=np.gradient(u2[it,ipres,:,iy],distlat)
        dvdy2=np.gradient(v2[it,ipres,ix,:],distlon)
    
    dh1a[it,:]=htad1[it,:]-adh1[it,:]   #sign htad is negative, adh pos due to down=pos
    dh1b[it,:]=dt1[it,:]+htad1[it,:]-adh1[it,:]
    dh2a[it,:]=htad2[it,:]-adh2[it,:]   #sign htad is negative, adh pos due to down=pos
    dh2b[it,:]=dt2[it,:]+htad2[it,:]-adh2[it,:]

    #Integrate over time
    if it>0:
        for ipres in range(0,npres):
            temp_adh1[it,ipres]=temp_adh1[it-1,ipres]+(adh1[it,ipres]+adh1[it-1,ipres])*1800.
            temp_htad1[it,ipres]=temp_htad1[it-1,ipres]+(htad1[it,ipres]+htad1[it-1,ipres])*1800.
            temp_dh1a[it,ipres]=temp_dh1a[it-1,ipres]+(dh1a[it,ipres]+dh1a[it-1,ipres])*1800.
            temp_dh1b[it,ipres]=temp_dh1b[it-1,ipres]+(dh1b[it,ipres]+dh1b[it-1,ipres])*1800.
                     
            temp_adh2[it,ipres]=temp_adh2[it-1,ipres]+(adh2[it,ipres]+adh2[it-1,ipres])*1800.
            temp_htad2[it,ipres]=temp_htad2[it-1,ipres]+(htad2[it,ipres]+htad2[it-1,ipres])*1800.
            temp_dh2a[it,ipres]=temp_dh2a[it-1,ipres]+(dh2a[it,ipres]+dh2a[it-1,ipres])*1800.
            temp_dh2b[it,ipres]=temp_dh2b[it-1,ipres]+(dh2b[it,ipres]+dh2b[it-1,ipres])*1800.

#Take the means over the entire three months                     
mean_dt1=np.mean(dt1,axis=0)    #per sec
mean_adh1=np.mean(adh1,axis=0)
mean_htad1=np.mean(htad1,axis=0)
mean_dh1a=np.mean(dh1a,axis=0)
mean_dh1b=np.mean(dh1b,axis=0)

mean_temp_adh1=np.mean(temp_adh1,axis=0)
mean_temp_htad1=np.mean(temp_htad1,axis=0)
mean_temp_dh1a=np.mean(temp_dh1a,axis=0)
mean_temp_dh1b=np.mean(temp_dh1b,axis=0)
                
mean_dt2=np.mean(dt2,axis=0)    #per sec
mean_adh2=np.mean(adh2,axis=0)

mean_temp_adh2=np.mean(temp_adh2,axis=0)
mean_temp_htad2=np.mean(temp_htad2,axis=0)
mean_temp_dh2a=np.mean(temp_dh2a,axis=0)
mean_temp_dh2b=np.mean(temp_dh2b,axis=0)

#Plot the means of the integrated temperatures
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'font.family': 'sans-serif'})
plt.rcParams.update({'font.weight': 'light'})

fig, ax = plt.subplots(1,1,sharex=True,figsize=(8,6))
plt.plot((mean_t2-mean_t1)*100.,hgt,c='grey',linestyle='--',linewidth=2.5)
plt.plot((mean_temp_adh2-mean_temp_adh1),hgt,c='r',linewidth=2.5)
plt.plot(-(mean_temp_htad2-mean_temp_htad1),hgt,c='g',linewidth=2.5)
plt.plot((mean_temp_dh2a-mean_temp_dh1a),hgt,c='b',linewidth=2.5)
plt.plot([0,0],[0,1000],c='k',linewidth=0.5)
#plt.plot((mean_temp_dh2b-mean_temp_dh1b),hgt,c='b',linestyle='--')
plt.plot((dt0)*100,hgt,c='orange',linewidth=2.5)
plt.ylim([0,np.max(hgt)])
ax.set_yticks([hgt[7],hgt[5],hgt[3],hgt[0]])
ax.set_yticklabels(preslabel)
ax.annotate('x100',xy=(0.08,0.8),xycoords='axes fraction',size=20,color='orange')
ax.annotate('x100',xy=(0.85,0.6),xycoords='axes fraction',size=20,color='grey')
plt.xlabel('Temperature difference (K)')
plt.ylabel('Pressure (hPa)')
plt.savefig('dekoketal_thermo.pdf',bbox_inches='tight') 