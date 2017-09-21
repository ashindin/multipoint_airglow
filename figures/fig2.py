
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


# In[5]:

bsc_fn="../astrometric_calibration/bsc5.dat"
cl_fn="../astrometric_calibration/ConstellationLines.dat"


# In[6]:

fid=open(bsc_fn,'r')
bsc_lines=fid.readlines()
fid.close()


# In[7]:

fid=open(cl_fn,'r')
cl_lines=fid.readlines()
cl_lines=cl_lines[7::]
cl_lines=[cl[4:-1] for cl in cl_lines]
fid.close()
xx_list_sbig=[]
yy_list_sbig=[]
fname="../data/160829/sbig/red1_378.FIT"
save_fname="../astrometric_calibration/sbig_160829_solve.pars"
lat_deg=56.1501667
lon_deg=46.1050833
hei_m=183.
err,az0,alt0,a,b,c,d=get_solve_pars(save_fname)
sbig_site=EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg, height=hei_m*u.m)
date_obs=sbig_get_date_obs(fname)

i=0
# for cl in cl_lines[46:47]:
for cl in cl_lines:
    temp=(cl.split()[1::])
#     print(temp)
    st_list=[int(t) for t in temp]
#     print(st_list)
    
    st_sao_list=[]
    for st in st_list:
#         print(st)
        try:
            st_sao=int(bsc_lines[st-1][31:37])
            st_sao_list.append(st_sao)
        except:
            continue
    
    ra_list=[]
    dec_list=[]
    
    for sao in st_sao_list:
#         print(sao)
        try:
            temp=np.where(sao_nums==sao)[0][0]
        except:
            continue
        ra_list.append(ra[temp])
        dec_list.append(dec[temp])
#         print(sao,ra[temp],dec[temp])
        
    
    SC=SkyCoord(ra_list, dec_list, frame='icrs', unit='deg');
    AA=SC.transform_to(AltAz(obstime=date_obs, location=sbig_site,temperature=20*u.deg_C,pressure=1013*u.hPa,
                                                   relative_humidity=0.5,obswl=630.0*u.nm))

    az_list=AA.az.rad
    alt_list=AA.alt.rad

#     print(alt_list*180/np.pi)
    
    x_list, y_list = tan_hor2pix(az_list,alt_list,az0,alt0,c,d)
    xx_list_sbig.append(x_list)
    yy_list_sbig.append(y_list)
    
#     print(i)
    i+=1
#     plt.plot(x_list,y_list,'r')


# In[9]:

fid=open(cl_fn,'r')
cl_lines=fid.readlines()
cl_lines=cl_lines[7::]
cl_lines=[cl[4:-1] for cl in cl_lines]
fid.close()
xx_list_keo=[]
yy_list_keo=[]
fname="../data/160829/keo/Capture_442.fit"
save_fname="../astrometric_calibration/keo_160829_solve_manual.pars"
lat_deg=55.9305361
lon_deg=48.7444861
hei_m=91.
err,az0,alt0,a,b,c,d=get_solve_pars(save_fname)
keo_site=EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg, height=hei_m*u.m)
date_obs=keo_get_date_obs(fname)

i=0
# for cl in cl_lines[46:47]:
for cl in cl_lines:
    temp=(cl.split()[1::])
#     print(temp)
    st_list=[int(t) for t in temp]
#     print(st_list)
    
    st_sao_list=[]
    for st in st_list:
#         print(st)
        try:
            st_sao=int(bsc_lines[st-1][31:37])
            st_sao_list.append(st_sao)
        except:
            continue
    
    ra_list=[]
    dec_list=[]
    
    for sao in st_sao_list:
#         print(sao)
        try:
            temp=np.where(sao_nums==sao)[0][0]
        except:
            continue
        ra_list.append(ra[temp])
        dec_list.append(dec[temp])
#         print(sao,ra[temp],dec[temp])
        
    
    SC=SkyCoord(ra_list, dec_list, frame='icrs', unit='deg');
    AA=SC.transform_to(AltAz(obstime=date_obs, location=keo_site,temperature=20*u.deg_C,pressure=1013*u.hPa,
                                                   relative_humidity=0.5,obswl=630.0*u.nm))

    az_list=AA.az.rad
    alt_list=AA.alt.rad

#     print(alt_list*180/np.pi)
    
    x_list, y_list = arc_hor2pix(az_list,alt_list,az0,alt0,c,d)
    xx_list_keo.append(x_list)
    yy_list_keo.append(y_list)
    
#     print(i)
    i+=1
#     plt.plot(x_list,y_list,'r')


# In[ ]:

def figure2(fit_fn1,fit_fn2,save_fname1,save_fname2, xx_list_sbig, yy_list_sbig, xx_list_keo, yy_list_keo, lat_deg_sbig=56.1501667,lon_deg_sbig=46.1050833,hei_m_sbig=183.,lat_deg_keo=55.9305361,lon_deg_keo=48.7444861,hei_m_keo=91.):

    
    err_sbig,az0_sbig,alt0_sbig,a_sbig,b_sbig,c_sbig,d_sbig=get_solve_pars(save_fname1)
    err_keo,az0_keo,alt0_keo,a_keo,b_keo,c_keo,d_keo=get_solve_pars(save_fname2)

    
    fig=plt.figure(figsize=(6.69,3.4))
    

    ax1=plt.axes(position=[0.06, 0.1+0.12, 0.43, 0.64])   
    ax2=plt.axes(position=[0.55, 0.1, 0.43, 0.85])
         

    plt.sca(ax1)
    vshift=2000
    fname=fit_fn1

    sbig_site=EarthLocation(lat=lat_deg_sbig*u.deg, lon=lon_deg_sbig*u.deg, height=hei_m_sbig*u.m)
    date_obs=sbig_get_date_obs(fname)

    md_fname='../spectrophotometric_calibration/sbig_160829_masterdark.fit'
    mf_fname='../spectrophotometric_calibration/sbig_master.flat'
    img=fits.getdata(fname).astype('float')
    masterdark=fits.getdata(md_fname).astype('float')
    masterflat=fits.getdata(mf_fname).astype('float')
    img=(img-masterdark)/masterflat

    vrange=1.5*np.std(img)
    pcm1=plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-vrange+vshift, vmax=np.median(img)+vrange+vshift)
    pcm1.set_edgecolor('face')
    plt.axis('equal')
    ax1.set_xlim((0,372))
    ax1.set_ylim((281,0)) 
    for m in [19,20,47]:
        x_list=xx_list_sbig[m]
        y_list=yy_list_sbig[m]
        plt.plot(x_list,y_list,'y',lw=1.5)

    az_axe=np.linspace(0,2*np.pi,100)
    for alt_deg in range(80,92,2):            
        alt_axe=np.pi/180*alt_deg*np.ones_like(az_axe)    
        x_axe, y_axe = tan_hor2pix(az_axe,alt_axe,az0_sbig,alt0_sbig,c_sbig,d_sbig)
        plt.plot(x_axe,y_axe,'g',lw=0.5)

    alt_axe=np.linspace(75*np.pi/180,88*np.pi/180,100)
    for az_deg in range(0,360,10):            
        az_axe=np.pi/180*az_deg*np.ones_like(alt_axe)    
        x_axe, y_axe = tan_hor2pix(az_axe,alt_axe,az0_sbig,alt0_sbig,c_sbig,d_sbig)
        if az_deg==180:
            plt.plot(x_axe,y_axe,'r',lw=0.5)
        else:
            plt.plot(x_axe,y_axe,'g',lw=0.5)

    plt.sca(ax2)
    vshift=300
    fname=fit_fn2
    keo_site=EarthLocation(lat=lat_deg_keo*u.deg, lon=lon_deg_keo*u.deg, height=hei_m_keo*u.m)
    date_obs=keo_get_date_obs(fname)
    md_fname='../spectrophotometric_calibration/keo_160829_masterdark.fit'        
    img=fits.getdata(fname).astype('float')
    masterdark=fits.getdata(md_fname).astype('float')        
    img=(img-masterdark)
    vrange=1.5*np.std(img)
    pcm2=plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-vrange+vshift, vmax=np.median(img)+vrange+vshift)
    pcm2.set_edgecolor('face')
    plt.axis('equal')

#         for m in [19,20,47]:
    for m in range(len(xx_list_keo)):
        x_list=xx_list_keo[m]
        y_list=yy_list_keo[m]
        plt.plot(x_list,y_list,'y',lw=0.7)

    az_axe=np.linspace(0,2*np.pi,100)
    for alt_deg in range(0,90,10):            
        alt_axe=np.pi/180*alt_deg*np.ones_like(az_axe)    
        x_axe, y_axe = arc_hor2pix(az_axe,alt_axe,az0_keo,alt0_keo,c_keo,d_keo)
        plt.plot(x_axe,y_axe,'g',lw=0.5)

    alt_axe=np.linspace(0*np.pi/180,80*np.pi/180,100)
    for az_deg in range(0,360,10):            
        az_axe=np.pi/180*az_deg*np.ones_like(alt_axe)    
        x_axe, y_axe = arc_hor2pix(az_axe,alt_axe,az0_keo,alt0_keo,c_keo,d_keo)
        if az_deg==180:
            plt.plot(x_axe,y_axe,'r',lw=0.5)
        else:
            plt.plot(x_axe,y_axe,'g',lw=0.5)

    ax2.set_xlim((511,0))
    ax2.set_ylim((511,0))
    
    
    plt.savefig('fig2.png',dpi=300)
    plt.savefig('fig2.eps',dpi=300)
    plt.savefig('fig2.pdf',dpi=300)
    plt.close()
#     plt.show()

fit_fn1="../data/160829/sbig/red1_378.FIT"
fit_fn2="../data/160829/keo/Capture_442.fit"
save_fname1="../astrometric_calibration/sbig_160829_solve.pars"
save_fname2="../astrometric_calibration/keo_160829_solve_manual.pars"
figure2(fit_fn1,fit_fn2,save_fname1,save_fname2,xx_list_sbig, yy_list_sbig, xx_list_keo, yy_list_keo)

