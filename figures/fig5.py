
# coding: utf-8

# In[1]:

import datetime
import os
import sys, inspect
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# In[2]:

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../astrometric_calibration")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from tan_module import *
from arc_module import *


# In[3]:

from astropy.utils import iers
iers.conf.iers_auto_url="http://toshi.nofs.navy.mil/ser7/finals2000A.all"
from astropy.time import Time
t = Time('2016:001')
t.ut1  


# In[4]:

def sbig_get_date_obs(filename):
    fid_fit=fits.open(filename);
    date_obs_str=fid_fit[0].header["DATE-OBS"]  
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S.%f")
def keo_get_date_obs(filename):
    fid_fit=fits.open(filename);
    date_obs_str=fid_fit[0].header["DATE-OBS"]  
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S")

def get_icrs_coord_of_bright_stars(filename='cat_icrs_sao_bright_35364.dat'):
    fid=open(filename,'r')
    Lines=fid.readlines()
    fid.close()
    num_stars=len(Lines)
    sao_nums=np.zeros(num_stars,dtype=np.int)
    ra=np.zeros(num_stars)
    dec=np.zeros(num_stars)
    mag=np.zeros(num_stars)
    for i in range(num_stars):
        par_list=Lines[i][4:-1].split(" ")
        sao_nums[i]=par_list[0]
        ra[i]=par_list[1]
        dec[i]=par_list[2]
        mag[i]=par_list[3]
    return sao_nums, ra, dec, mag
sao_nums, ra, dec, mag = get_icrs_coord_of_bright_stars(filename='../astrometric_calibration/cat_icrs_sao_bright_35364.dat')

def get_solve_pars(save_fname):
    a=np.zeros(3)
    b=np.zeros(3)
    c=np.zeros(3)
    d=np.zeros(3)
    fid=open(save_fname,'r')
    lines=fid.readlines()
    pars_list=lines[1].split(' ')
    pars=[float(par_str) for par_str in pars_list]
    err,az0,alt0,a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2],d[0],d[1],d[2]=pars
    fid.close()
    return err,az0,alt0,a,b,c,d


# In[59]:

def figure2(fit_fn1,fit_fn2):

    
#     err_sbig,az0_sbig,alt0_sbig,a_sbig,b_sbig,c_sbig,d_sbig=get_solve_pars(save_fname1)
#     err_keo,az0_keo,alt0_keo,a_keo,b_keo,c_keo,d_keo=get_solve_pars(save_fname2)

    
    fig=plt.figure(figsize=(6.69,3.4))
    

    ax1=plt.axes(position=[0.055, 0.1, 0.4, 0.79])   
    ax2=plt.axes(position=[0.55-0.035, 0.1, 0.4, 0.79])


    plt.sca(ax1)
    
    fname=fit_fn1

#     sbig_site=EarthLocation(lat=lat_deg_sbig*u.deg, lon=lon_deg_sbig*u.deg, height=hei_m_sbig*u.m)
    date_obs=keo_get_date_obs(fname)

    md_fname='../spectrophotometric_calibration/keo_140826_masterdark.fit'        
    img=fits.getdata(fname).astype('float')
    masterdark=fits.getdata(md_fname).astype('float')        
    img=(img-masterdark)

    vshift=0
    vrange=2*np.std(img[334:434,186:286])
    img_med=np.median(img[334:434,186:286])
    pcm1=plt.pcolormesh(img, cmap="gray", vmin=img_med-vrange+vshift, vmax=img_med+vrange+vshift)
    pcm1.set_edgecolor('face')
    plt.axis('equal')
    ax1.set_xlim((511,0))
    ax1.set_ylim((511,0)) 

    plt.plot([200-14,300-14],[325+9,325+9],'r',lw=2)
    plt.plot([200-14,300-14],[425+9,425+9],'r',lw=2)    
    plt.plot([200-14,200-14],[325+9,425+9],'r',lw=2)
    plt.plot([300-14,300-14],[325+9,425+9],'r',lw=2)
    
    plt.sca(ax2)    
    fname=fit_fn2
#     keo_site=EarthLocation(lat=lat_deg_keo*u.deg, lon=lon_deg_keo*u.deg, height=hei_m_keo*u.m)
#     date_obs=keo_get_date_obs(fname)    
    img=fits.getdata(fname).astype('float')
#     masterdark=fits.getdata(md_fname).astype('float')        
#     img=(img-masterdark)
#     vrange=1.5*np.std(img)
    pcm2=plt.pcolormesh(img, cmap="jet", vmin=0, vmax=20)
    pcm2.set_edgecolor('face')
    plt.axis('equal')


    ax2.set_xlim((300-14,200-14))
    ax2.set_ylim((425+9,325+9)) 
    
#     ax2.set_yticklabels((""))
    
    ax2_cb=plt.axes(position=[0.55+0.43+0.008-0.065, 0.1, 0.015, 0.79])
    plt.colorbar(pcm2,ax2_cb)  
    
#     plt.plot([200-14,300-14],[325+9,325+9],'r',lw=2)
#     plt.plot([200-14,300-14],[425+9,425+9],'r',lw=2)    
#     plt.plot([200-14,200-14],[325+9,425+9],'r',lw=2)
#     plt.plot([300-14,300-14],[325+9,425+9],'r',lw=2)
    
    
    plt.savefig('fig5.png',dpi=300)
    plt.savefig('fig5.eps')
    plt.savefig('fig5.pdf')
#     plt.close()
    #plt.show()

fit_fn1="../data/140826/keo/20140826_3214.fit"
fit_fn2="../image_processing/keo140826_glowfit/20140826_3214.fit"

figure2(fit_fn1,fit_fn2)

