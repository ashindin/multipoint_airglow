
# coding: utf-8

# In[1]:

import datetime
import os
import sys
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from arc_module import *


# In[5]:

fit_fname="../data/140824/keo/20140824_2400.fit"
save_fname="keo_140824_solve.pars"
#save_fname="keo_140824_solve_manual.pars"
X0=230
Y0=380
# X0=255
# Y0=255


# In[6]:

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

sao_nums, ra, dec, mag = get_icrs_coord_of_bright_stars()
num_stars=500
sao_nums_filt = sao_nums[0:num_stars]
ra_filt= ra[0:num_stars]
dec_filt= dec[0:num_stars]
mag_filt= mag[0:num_stars]

def get_solve_pars(save_fname):
    a=np.zeros(3)
    b=np.zeros(3)
    c=np.zeros(3)
    d=np.zeros(3)
    fid=open(save_fname,'r')
    lines=fid.readlines()
    pars_list=lines[1].split(' ')
    pars=[float(par_str) for par_str in pars_list]
    az0,alt0,a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2],d[0],d[1],d[2]=pars
    fid.close()
    return az0,alt0,a,b,c,d


# In[7]:

def keo_manual_solving(fit_fname, save_fname, lat_deg=55.9305361,lon_deg=48.7444861,hei_m=91.):
    
    az0,alt0,a,b,c,d=get_solve_pars(save_fname)
        
    ms=5
    mew=1
    
    fig=plt.figure(figsize=(12.8,7.2))
    fig.set_size_inches(12.8, 7.2)
    
    dx=0.05    
    ax=plt.axes(position=[0.035000000000000003+dx/2, 0.061764705882352888, 0.50624999999999998, 0.89999999999999991])   
    plt.axis('equal')
    ax1=plt.axes(position=[0.58205882352941185+dx, 0.71052941176470596, 0.140625, 0.25])
    ax2=plt.axes(position=[0.7720588235294118+dx, 0.71052941176470596, 0.140625, 0.25])
    ax3=plt.axes(position=[0.58205882352941185+dx, 0.39552941176470596, 0.140625, 0.24999999999999994])
    ax4=plt.axes(position=[0.7720588235294118+dx, 0.39552941176470596, 0.140625, 0.24999999999999994])
    ax5=plt.axes(position=[0.58205882352941185+dx, 0.061764705882352888, 0.140625, 0.25])
    ax6=plt.axes(position=[0.7720588235294118+dx, 0.061764705882352888, 0.140625, 0.25])
                 
    fname=fit_fname
    keo_site=EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg, height=hei_m*u.m)
    date_obs=keo_get_date_obs(fname)

    SC=SkyCoord(ra_filt, dec_filt, frame='icrs', unit='deg');

    AA=SC.transform_to(AltAz(obstime=date_obs, location=keo_site,temperature=20*u.deg_C,pressure=1013*u.hPa,
                                               relative_humidity=0.5,obswl=630.0*u.nm))

    AZ=AA.az.rad
    ALT=AA.alt.rad

    filt_mask=np.where(ALT>30.0*np.pi/180)
    ALT=ALT[filt_mask]
    AZ=AZ[filt_mask]
    sao_nums_filt2=sao_nums_filt[filt_mask]

    X,Y = arc_hor2pix(AZ,ALT,az0,alt0,c,d)

    R=np.sqrt((X-X0)**2+(Y-Y0)**2)
    I=np.argsort(R)
    X=X[I]
    Y=Y[I]
    sao_nums_filt2=sao_nums_filt2[I]
    AZ=AZ[I]
    ALT=ALT[I]

    hdulist = fits.open(fname,ignore_missing_end=True)
    img=hdulist[0].data
    hdulist.close()
    plt.sca(ax)
    plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-300, vmax=np.median(img)+300)
    plt.plot(X[0],Y[0],marker="o", lw=0.,mec="g", mfc='none',mew=mew, ms=ms)
    plt.text(X[0], Y[0], str(sao_nums_filt2[0]))
    for j in range(1,6):
        plt.plot(X[j],Y[j],marker="o", lw=0.,mec="b",mew=mew, ms=ms, mfc='none')
        plt.text(X[j], Y[j], str(sao_nums_filt2[j]))
    plt.plot(X[5::],Y[5::],marker="o", lw=0.,mew=mew,mec="r", mfc='none',ms=ms)
    for j in range(6,len(sao_nums_filt2)):
        plt.text(X[j], Y[j], str(sao_nums_filt2[j]),fontsize="smaller")
    ax.set_xlim((511,1))
    ax.set_ylim((511,1))
    plt.title(fname.split('/')[-1] + " " + "{0:0>2}".format(date_obs.hour) + ":" + "{0:0>2}".format(date_obs.minute) + ":" + "{0:0>2}".format(date_obs.second))

    plt.sca(ax1)
    plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-300, vmax=np.median(img)+300)
    plt.plot(X[0],Y[0],color="g",marker="o", lw=0.,mec="g", mfc='none',mew=mew, ms=ms)
    plt.title("SAO " + str(sao_nums_filt2[0]))
    ax1.set_xlim((X[0]+20,X[0]-20))
    ax1.set_ylim((Y[0]+20,Y[0]-20))

    plt.sca(ax2)
    plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-300, vmax=np.median(img)+300)
    plt.plot(X[1],Y[1],color="b",marker="o", lw=0.,mec="b", mfc='none',mew=mew, ms=ms)
    plt.title("SAO " + str(sao_nums_filt2[1]))
    ax2.set_xlim((X[1]+20,X[1]-20))
    ax2.set_ylim((Y[1]+20,Y[1]-20))

    plt.sca(ax3)
    plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-300, vmax=np.median(img)+300)
    plt.plot(X[2],Y[2],color="b",marker="o", lw=0.,mec="b", mfc='none',mew=mew, ms=ms)
    plt.title("SAO " + str(sao_nums_filt2[2]))
    ax3.set_xlim((X[2]+20,X[2]-20))
    ax3.set_ylim((Y[2]+20,Y[2]-20))

    plt.sca(ax4)
    plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-300, vmax=np.median(img)+300)
    plt.plot(X[3],Y[3],color="b",marker="o", lw=0.,mec="b", mfc='none',mew=mew, ms=ms)
    plt.title("SAO " + str(sao_nums_filt2[3]))
    ax4.set_xlim((X[3]+20,X[3]-20))
    ax4.set_ylim((Y[3]+20,Y[3]-20))

    plt.sca(ax5)
    plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-300, vmax=np.median(img)+300)
    plt.plot(X[4],Y[4],color="b",marker="o", lw=0.,mec="b", mfc='none',mew=mew, ms=ms)
    plt.title("SAO " + str(sao_nums_filt2[4]))
    ax5.set_xlim((X[4]+20,X[4]-20))
    ax5.set_ylim((Y[4]+20,Y[4]-20))

    plt.sca(ax6)
    plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-300, vmax=np.median(img)+300)
    plt.plot(X[5],Y[5],color="b",marker="o", lw=0.,mec="b", mfc='none',mew=mew, ms=ms)
    plt.title("SAO " + str(sao_nums_filt2[5]))
    ax6.set_xlim((X[5]+20,X[5]-20))
    ax6.set_ylim((Y[5]+20,Y[5]-20))
    plt.show()
keo_manual_solving(fit_fname, save_fname)

