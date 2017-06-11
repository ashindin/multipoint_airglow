
# coding: utf-8

# In[1]:

import datetime
import os, sys, inspect
import numpy as np
from astropy.io import fits
import scipy.signal as ss

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../astrometric_calibration")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from tan_module import *
from matplotlib import dates


# In[2]:

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)


# In[3]:

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz


# In[4]:

def s1c_get_date_obs(filename,ut_shift=-3):
    fdt=open(filename,'r',errors = 'ignore')
    fdt.seek(8*80)
    test=fdt.read(39)
    fdt.close()
    test=test[15::]
    ms=test[16:20]
    if ms=='1000':
        test=test[0:16]+'000 2014'
#         test[16]='0'
        obs_time=datetime.datetime.strptime(test,'%b %d %H:%M:%S.%f %Y')+datetime.timedelta(hours=ut_shift)+datetime.timedelta(seconds=1)
    else:
        obs_time=datetime.datetime.strptime(test,'%b %d %H:%M:%S.%f %Y')+datetime.timedelta(hours=ut_shift)
    return obs_time

def s1c_get_exp_sec(filename):
    fdt=open(filename,'r',errors = 'ignore')
    fdt.seek(7*80)
    temp=fdt.read(80)
    exp_sec=float(temp.split("=")[-1].split(" ")[1])
    fdt.close()
    return exp_sec

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

def get_spcal_day_coefs(spcal_day_fname):
    fid = open(spcal_day_fname,'r')
    lines=fid.readlines()
    fid.close()
    lines_list=lines[1].split(" ")
    spcal_coef=float(lines_list[0])
    spcal_std=float(lines_list[1])
    return spcal_coef, spcal_std


# In[18]:

s1c_fit_path="../data/140824/s1c"
base_frames_fname="s1c_140824_base.frames"
solve_pars_fname="../astrometric_calibration/s1c_140824_solve.pars"
s1c_spcal_day_fname = "../spectrophotometric_calibration/s1c_140824_day.spcal"
masterdark_fname="../spectrophotometric_calibration/s1c_140824_masterdark.fit"
masterflat_fname="../spectrophotometric_calibration/s1c_master.flat"
lat_cam_deg=56.1501667; lat_cam=lat_cam_deg*np.pi/180;
lon_cam_deg=46.1050833; lon_cam=lon_cam_deg*np.pi/180;
hei_cam=183.;
CAM_site=EarthLocation(lat=lat_cam_deg*u.deg, lon=lon_cam_deg*u.deg, height=hei_cam*u.m)
avr_width1=11
avr_width2=59
interp_deg=2 # three points


# In[6]:

s1c_spcal_day_coef, s1c_spcal_std=get_spcal_day_coefs(s1c_spcal_day_fname)   


# In[7]:

spath="./s1c140824_glowfit/"
if not os.path.exists(spath):
    os.makedirs(spath)


# In[8]:

hdulist = fits.open(masterdark_fname,ignore_missing_end=True)
masterdark=hdulist[0].data
hdulist.close()
hdulist = fits.open(masterflat_fname,ignore_missing_end=True)
masterflat=hdulist[0].data
hdulist.close()


# In[9]:

masterdark.shape


# In[10]:

fid=open(base_frames_fname,'r')
lines=fid.readlines()
fid.close()
base_frames=[]
for i in range(len(lines)):
    if lines[i]!='\n':               
        base_frames.append(lines[i][0:-1])
s1c_fit_filenames=sorted([s1c_fit_path+'/'+fn for fn in next(os.walk(s1c_fit_path))[2]])
base_frames_fullnames=[s1c_fit_path + "/" + bf_name for bf_name in base_frames]


# In[11]:

fr_inds=list(range(len(s1c_fit_filenames)))
bf_inds=[]
bf_dates=[]
bf_x_dates=[]
for bfn in base_frames_fullnames:
    for i in fr_inds:
        if bfn==s1c_fit_filenames[i]:
            bf_inds.append(i)
    bf_dates.append(s1c_get_date_obs(bfn,-4)+datetime.timedelta(seconds=s1c_get_exp_sec(bfn)/2))
    bf_x_dates.append(dates.date2num(bf_dates[-1]))
bf_inds_local=list(range(len(bf_inds)))
# bf_dates


# In[12]:

BF_imgs=np.zeros((masterdark.shape[0],masterdark.shape[1],len(base_frames_fullnames)))
for i in range(len(base_frames_fullnames)):
#     sys.stdout.write('\r')
#     sys.stdout.write("Add base frame "+str(i+1)+"/"+str(len(base_frames_fullnames)))
#     sys.stdout.flush()
    bf_fname=base_frames_fullnames[i]
    hdulist = fits.open(bf_fname,ignore_missing_end=True)
    BF_imgs[:,:,i]=hdulist[0].data.astype('float')
    hdulist.close()
    BF_imgs[:,:,i]=(BF_imgs[:,:,i]-masterdark.astype('float'))/masterflat.astype('float')
    BF_imgs[:,:,i]=ss.medfilt(BF_imgs[:,:,i],kernel_size=avr_width1)


# In[13]:

YPIX, XPIX = np.mgrid[1:BF_imgs[:,:0].shape[0]+1, 1:BF_imgs[:,:,0].shape[1]+1]


# In[14]:

def shift_img(img_base, base_time, obs_time):
    AZ=np.zeros_like(img_base)
    ALT=np.zeros_like(img_base)
    
    XPIX2=np.copy(XPIX)
    YPIX2=np.copy(YPIX)
    
    AZ,ALT=tan_pix2hor(XPIX,YPIX,az0,alt0,a,b)
    C=SkyCoord(alt = ALT*u.rad, az = AZ*u.rad, obstime = base_time, frame = 'altaz', location = CAM_site, temperature=15*u.deg_C,pressure=1013*u.hPa, relative_humidity=0.63,obswl=630.0*u.nm)
    ALTAZ_new=C.transform_to(AltAz(obstime = obs_time, location = CAM_site, temperature=15*u.deg_C,pressure=1013*u.hPa, relative_humidity=0.63,obswl=630.0*u.nm))
    XPIX2, YPIX2=tan_hor2pix(ALTAZ_new.az.rad,ALTAZ_new.alt.rad,az0,alt0,c,d)

    med_value=np.median(img_base)
    img_mod=np.zeros_like(img_base)
    for i in range(np.size(img_base)):
        x=XPIX2.flat[i]-np.trunc(XPIX2.flat[i]);
        y=YPIX2.flat[i]-np.trunc(YPIX2.flat[i]);
        x0=int(np.trunc(XPIX2.flat[i]));
        y0=int(np.trunc(YPIX2.flat[i]));
        if (x0>=1) and (y0>=1) and (x0<img_base.shape[1]) and (y0<img_base.shape[0]):
            img_mod[y0,x0]=img_mod[y0,x0]+(1-x)*(1-y)*img_base.flat[i];
        if (x0+1>=1) and (y0>=1) and (x0+1<img_base.shape[1]) and (y0<img_base.shape[0]):
            img_mod[y0,x0+1]=img_mod[y0,x0+1]+(x)*(1-y)*img_base.flat[i];
        if (x0>=1) and (y0+1>=1) and (x0<img_base.shape[1]) and (y0+1<img_base.shape[0]):
            img_mod[y0+1,x0]=img_mod[y0+1,x0]+(1-x)*(y)*img_base.flat[i];
        if (x0+1>=1) and (y0+1>=1) and (x0+1<img_base.shape[1]) and (y0+1<img_base.shape[0]):
            img_mod[y0+1,x0+1]=img_mod[y0+1,x0+1]+(x)*(y)*img_base.flat[i];
    img_mod.flat[np.where(img_mod.flat == 0.)[0]]=med_value    
#     img_mod = si.griddata((YPIX2.flat,XPIX2.flat), img_base.flat, (XPIX, YPIX), fill_value=med_value, method='linear')
    return img_mod;


# In[15]:

az0,alt0,a,b,c,d=get_solve_pars(solve_pars_fname)


# In[48]:

left_bf_ind=0
right_bf_ind=0
for i in range(bf_inds[0],bf_inds[-1]):
# for i in range(36,37):
    sys.stdout.write('\r')
    sys.stdout.write("Processing frame "+str(i+1)+"/"+str(len(range(bf_inds[0],bf_inds[-1]))))
    sys.stdout.flush()
    bf_locals=[]
    bfd_locals=[]
    bfl_locals=[]
    bfdx_locals=[]
    fn=s1c_fit_filenames[i]
    f_exp=s1c_get_exp_sec(fn)
    f_date=s1c_get_date_obs(fn,-4)+datetime.timedelta(seconds=f_exp/2)
    f_x_date=dates.date2num(f_date)
    if i in bf_inds: continue;
    for j in bf_inds_local:
        if f_x_date>bf_x_dates[j]:
            bf_locals
            left_bf_ind=j
            right_bf_ind=j+1
    bf_locals.append(bf_inds[left_bf_ind])
    bf_locals.append(bf_inds[right_bf_ind])
    bfd_locals.append(bf_dates[left_bf_ind])
    bfd_locals.append(bf_dates[right_bf_ind])
    bfdx_locals.append(-f_x_date+bf_x_dates[left_bf_ind])
    bfdx_locals.append(-f_x_date+bf_x_dates[right_bf_ind])
    bfl_locals.append(left_bf_ind)
    bfl_locals.append(right_bf_ind)
    if left_bf_ind < bf_inds_local[-1]/2:
        for j in range(interp_deg-1):
            bf_locals.append(bf_inds[right_bf_ind+j+1])
            bfd_locals.append(bf_dates[right_bf_ind+j+1])
            bfdx_locals.append(-f_x_date+bf_x_dates[right_bf_ind+j+1])
            bfl_locals.append(right_bf_ind+j+1)
    else:
        for j in range(interp_deg-1):
            bf_locals.insert(0,bf_inds[left_bf_ind-j-1])        
            bfd_locals.insert(0,bf_dates[left_bf_ind-j-1])
            bfdx_locals.insert(0,-f_x_date+bf_x_dates[left_bf_ind-j-1])
            bfl_locals.insert(0,left_bf_ind-j-1)
#     print(i," - ", bf_locals,bfl_locals, bfdx_locals)

    y=[]
    for j in range(len(bfl_locals)):
        y.append(shift_img(BF_imgs[:,:,bfl_locals[j]],bfd_locals[j],f_date))
    
    det=bfdx_locals[0]**2*bfdx_locals[1]+bfdx_locals[2]**2*bfdx_locals[0]+bfdx_locals[1]**2*bfdx_locals[2]-bfdx_locals[2]**2*bfdx_locals[1]-bfdx_locals[1]**2*bfdx_locals[0]-bfdx_locals[0]**2*bfdx_locals[2]
    dark=(bfdx_locals[0]**2*bfdx_locals[1]*y[2]+bfdx_locals[2]**2*bfdx_locals[0]*y[1]+bfdx_locals[1]**2*bfdx_locals[2]*y[0]
        -bfdx_locals[2]**2*bfdx_locals[1]*y[0]-bfdx_locals[1]**2*bfdx_locals[0]*y[2]-bfdx_locals[0]**2*bfdx_locals[2]*y[1])/det
    
    hdulist = fits.open(fn,ignore_missing_end=True)
    img=hdulist[0].data.astype('float')
    img=ss.medfilt((img-masterdark.astype('float'))/masterflat.astype('float') - dark,kernel_size=avr_width2) * s1c_spcal_day_coef
    hdulist[0].data=img
    hdulist[0].header['BITPIX']=-64
    hdulist[0].header['EXPTIME']=f_exp
    hdulist.writeto(spath+fn.split('/')[-1],overwrite=True,output_verify='fix')
    hdulist.close()  

