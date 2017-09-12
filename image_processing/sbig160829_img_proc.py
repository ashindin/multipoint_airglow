
# coding: utf-8

# In[1]:

import datetime
import os, sys, inspect
import numpy as np
from astropy.io import fits
import scipy.signal as ss
import scipy.interpolate as si

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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import dates


# In[5]:

def sbig_get_date_obs(filename):
    fid_fit=fits.open(filename);
    date_obs_str=fid_fit[0].header["DATE-OBS"]  
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S.%f")

def sbig_get_exp_sec(filename):
    fid_fit=fits.open(filename);
    exp_sec=fid_fit[0].header["EXPTIME"]  
    fid_fit.close()
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
    err,az0,alt0,a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2],d[0],d[1],d[2]=pars
    fid.close()
    return err,az0,alt0,a,b,c,d

def get_spcal_day_coefs(spcal_day_fname):
    fid = open(spcal_day_fname,'r')
    lines=fid.readlines()
    fid.close()
    lines_list=lines[1].split(" ")
    spcal_coef=float(lines_list[0])
    spcal_std=float(lines_list[1])
    return spcal_coef, spcal_std

def get_pumping_scheme_from_file(pumping_scheme_file):
    pumping_scheme=[]
    fid=open(pumping_scheme_file,'r')
    lines=fid.readlines()[1::]
    lines2=[]
    for line in lines:
        if line[-1]=='\n':
            lines2.append(line[:-1])
        else:
            lines2.append(line)
        line_list=lines2[-1].split(',')
        pumping_scheme.append([datetime.datetime.strptime(line_list[0], "%Y.%m.%dT%H:%M:%S"),float(line_list[1]),float(line_list[2]),float(line_list[3]),float(line_list[4]),float(line_list[5])])
    fid.close()
    return pumping_scheme

def make_clean_pumping_scheme(pumping_scheme):
    clean_pumping=[]
    date_axe=[]
    x_dates=[]
    for ps in pumping_scheme:
        st_time=ps[0]
        tau=ps[1]
        T=ps[2]
        num_p=ps[3]

        num_sec=int(T*num_p*60)
        date_axe_ps=[st_time+datetime.timedelta(seconds=i) for i in range(-1,num_sec)]
        x_dates_ps=dates.date2num(date_axe_ps)

        duty_frac=tau/T
        pump=np.zeros(num_sec+1)
        for i in range(-1,num_sec):
            min_float=i/60;
            period_float=i/(T*60)
            period_offset=period_float-int(period_float)
            if period_offset<=duty_frac and period_offset>0:
                pump[i]=ps[4]

        clean_pumping.append(pump)
        date_axe.append(date_axe_ps)
        x_dates.append(x_dates_ps)
    return date_axe, x_dates, clean_pumping

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
    #~ img_mod=np.zeros_like(img_base)
    #~ for i in range(np.size(img_base)):
        #~ x=XPIX2.flat[i]-np.trunc(XPIX2.flat[i]);
        #~ y=YPIX2.flat[i]-np.trunc(YPIX2.flat[i]);
        #~ x0=int(np.trunc(XPIX2.flat[i]));
        #~ y0=int(np.trunc(YPIX2.flat[i]));
        #~ if (x0>=1) and (y0>=1) and (x0<img_base.shape[1]) and (y0<img_base.shape[0]):
            #~ img_mod[y0,x0]=img_mod[y0,x0]+(1-x)*(1-y)*img_base.flat[i];
        #~ if (x0+1>=1) and (y0>=1) and (x0+1<img_base.shape[1]) and (y0<img_base.shape[0]):
            #~ img_mod[y0,x0+1]=img_mod[y0,x0+1]+(x)*(1-y)*img_base.flat[i];
        #~ if (x0>=1) and (y0+1>=1) and (x0<img_base.shape[1]) and (y0+1<img_base.shape[0]):
            #~ img_mod[y0+1,x0]=img_mod[y0+1,x0]+(1-x)*(y)*img_base.flat[i];
        #~ if (x0+1>=1) and (y0+1>=1) and (x0+1<img_base.shape[1]) and (y0+1<img_base.shape[0]):
            #~ img_mod[y0+1,x0+1]=img_mod[y0+1,x0+1]+(x)*(y)*img_base.flat[i];
    #~ img_mod.flat[np.where(img_mod.flat == 0.)[0]]=med_value
    img_mod = si.griddata((YPIX2.flat,XPIX2.flat), img_base.flat, (XPIX, YPIX), fill_value=med_value, method='linear')
    return img_mod;


# In[6]:

sbig_fit_path="../data/160829/sbig"
base_frames_fname="sbig_160829_base.frames"
solve_pars_fname="../astrometric_calibration/sbig_160829_solve.pars"
sbig_spcal_day_fname = "../spectrophotometric_calibration/sbig_160829_day.spcal"
sbig_spcal_fname = "../spectrophotometric_calibration/sbig_160829.spcal"
masterdark_fname="../spectrophotometric_calibration/sbig_160829_masterdark.fit"
masterflat_fname="../spectrophotometric_calibration/sbig_master.flat"
lat_cam_deg=56.1501667; lat_cam=lat_cam_deg*np.pi/180;
lon_cam_deg=46.1050833; lon_cam=lon_cam_deg*np.pi/180;
hei_cam=183.;
CAM_site=EarthLocation(lat=lat_cam_deg*u.deg, lon=lon_cam_deg*u.deg, height=hei_cam*u.m)
avr_width1=15
avr_width2=31
interp_deg=1 


# In[7]:

sbig_spcal_day_coef, sbig_spcal_std=get_spcal_day_coefs(sbig_spcal_day_fname)
sbig_spcal_day_coef, sbig_spcal_std


# In[8]:

fid=open(sbig_spcal_fname,'r')
lines=fid.readlines()
spcal_fnames=[]
spcal_coefs=[]
for i in range(1,len(lines)):
    line_list=lines[i][:-1].split(' ')
    spcal_fnames.append(line_list[0])
    spcal_coefs.append(float(line_list[1]))
fid.close()


# In[9]:

spath="./sbig160829_glowfit/"
if not os.path.exists(spath):
    os.makedirs(spath)


# In[10]:

hdulist = fits.open(masterdark_fname,ignore_missing_end=True)
masterdark=hdulist[0].data
hdulist.close()
hdulist = fits.open(masterflat_fname,ignore_missing_end=True)
masterflat=hdulist[0].data
hdulist.close()


# In[11]:

fid=open(base_frames_fname,'r')
lines=fid.readlines()
fid.close()
base_frames=[]
for i in range(len(lines)):
    if lines[i]!='\n':
        base_frames.append(lines[i][0:-1])
sbig_fit_filenames=sorted([sbig_fit_path+'/'+fn for fn in next(os.walk(sbig_fit_path))[2]])
base_frames_fullnames=[sbig_fit_path + "/" + bf_name for bf_name in base_frames]


# In[12]:

fr_inds=list(range(len(sbig_fit_filenames)))
bf_inds=[]
bf_dates=[]
bf_x_dates=[]
for bfn in base_frames_fullnames:
    for i in fr_inds:
        if bfn==sbig_fit_filenames[i]:
            bf_inds.append(i)
    bf_dates.append(sbig_get_date_obs(bfn)+datetime.timedelta(seconds=sbig_get_exp_sec(bfn)/2))
    bf_x_dates.append(dates.date2num(bf_dates[-1]))
bf_inds_local=list(range(len(bf_inds)))
# bf_dates


# In[13]:

# BF_imgs=np.zeros((masterdark.shape[0],masterdark.shape[1],len(base_frames_fullnames)))
BF_imgs=[]
for i in range(len(base_frames_fullnames)):
#     sys.stdout.write('\r')
#     sys.stdout.write("Add base frame "+str(i+1)+"/"+str(len(base_frames_fullnames)))
#     sys.stdout.flush()
    bf_fname=base_frames_fullnames[i]
    hdulist = fits.open(bf_fname,ignore_missing_end=True)
#     BF_imgs[:,:,i]=hdulist[0].data.astype('float')
    BF_imgs.append(hdulist[0].data.astype('float'))
    hdulist.close()
#     BF_imgs[:,:,i]=(BF_imgs[:,:,i]-masterdark.astype('float'))/masterflat.astype('float')
#     BF_imgs[:,:,i]=ss.medfilt(BF_imgs[:,:,i],kernel_size=avr_width1)
    BF_imgs[-1]=(BF_imgs[-1]-masterdark.astype('float'))/masterflat.astype('float')
    BF_imgs[-1]=ss.medfilt(BF_imgs[-1],kernel_size=avr_width1)


# In[14]:

# YPIX, XPIX = np.mgrid[1:BF_imgs[:,:,0].shape[0]+1, 1:BF_imgs[:,:,0].shape[1]+1]
YPIX, XPIX = np.mgrid[1:BF_imgs[0].shape[0]+1, 1:BF_imgs[0].shape[1]+1]


# In[15]:

err,az0,alt0,a,b,c,d=get_solve_pars(solve_pars_fname)


# In[16]:

pumping_scheme_file="pump160829.scheme"
prohibited_area_min=2.
pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)
pumping_scheme_list = [', '.join((ps[0].strftime('%Y-%m-%dT%H:%M:%S'),str(ps[1]) + ' min',str(ps[2])+' min',str(ps[3]),str(ps[4]),str(ps[5])+' kHz')) for ps in pumping_scheme]
date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)


# In[40]:

left_bf_ind=0
right_bf_ind=0
for i in range(bf_inds[0],bf_inds[-1]):
# for i in range(36,37):
    sys.stdout.write('\r')
    sys.stdout.write("Processing frame "+str(i+1)+"/"+str(len(range(bf_inds[-1]))))
    sys.stdout.flush()
    bf_locals=[]
    bfd_locals=[]
    bfl_locals=[]
    bfdx_locals=[]
    fn=sbig_fit_filenames[i]
    f_exp=sbig_get_exp_sec(fn)
    f_date_start=sbig_get_date_obs(fn)
    f_date=f_date_start+datetime.timedelta(seconds=f_exp/2)
    f_date_end=f_date_start+datetime.timedelta(seconds=f_exp)
    f_x_date_start=dates.date2num(f_date_start)
    f_x_date=dates.date2num(f_date)
    f_x_date_end=dates.date2num(f_date_end)
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

    fn_split=fn.split('/')[-1]
    sp_coef=0.
    sp_offset=0.
    for j in range(len(spcal_fnames)):
        if fn_split==spcal_fnames[j]:
            sp_coef=spcal_coefs[j]
            if sp_coef!=0.0:
                sp_offset=abs(sbig_spcal_day_coef-sp_coef)/sbig_spcal_std

    y=[]
    y_date=[]
    for j in range(len(bfl_locals)):
        y.append(shift_img(BF_imgs[bfl_locals[j]],bfd_locals[j],f_date))
        y_date.append(bfd_locals[j])
#     print (y_date)

    date_range=dates.date2num(y_date[1])-dates.date2num(y_date[0])
    dfraq2=(dates.date2num(f_date)-dates.date2num(y_date[0]))/date_range
    dfraq1=1-dfraq2
    #~ print(dfraq1,dfraq2)
    
#     det=bfdx_locals[0]**2*bfdx_locals[1]+bfdx_locals[2]**2*bfdx_locals[0]+bfdx_locals[1]**2*bfdx_locals[2]-bfdx_locals[2]**2*bfdx_locals[1]-bfdx_locals[1]**2*bfdx_locals[0]-bfdx_locals[0]**2*bfdx_locals[2]
#     dark=(bfdx_locals[0]**2*bfdx_locals[1]*y[2]+bfdx_locals[2]**2*bfdx_locals[0]*y[1]+bfdx_locals[1]**2*bfdx_locals[2]*y[0]
#         -bfdx_locals[2]**2*bfdx_locals[1]*y[0]-bfdx_locals[1]**2*bfdx_locals[0]*y[2]-bfdx_locals[0]**2*bfdx_locals[2]*y[1])/det
    dark=dfraq1*y[0]+dfraq2*y[1]

    f_date_iso=f_date_start.strftime('%Y-%m-%dT%H:%M:%S.%f')

    hdu_dark = fits.PrimaryHDU(dark*sbig_spcal_day_coef)
    hdu_dark.header['DATE-OBS']=f_date_iso
    hdu_dark.header['BITPIX']=-64
    hdu_dark.header['EXPTIME']=f_exp
    hdu_dark.header['PROJ-TYP']='TAN'
    hdu_dark.header['AMCAL-A0']=az0
    hdu_dark.header['AMCAL-H0']=alt0
    hdu_dark.header['AMCAL-A']=str(a)
    hdu_dark.header['AMCAL-B']=str(b)
    hdu_dark.header['AMCAL-C']=str(c)
    hdu_dark.header['AMCAL-D']=str(d)
    hdu_dark.header['MED-WID1']=avr_width1
    hdulist_dark = fits.HDUList([hdu_dark])
    hdulist_dark.writeto(spath+fn.split('/')[-1].split('.')[-2]+"_dark.fit",overwrite=True)

    hdulist = fits.open(fn,ignore_missing_end=True)
    img=hdulist[0].data.astype('float')
    #~ img=ss.medfilt((img-masterdark.astype('float'))/masterflat.astype('float') - dark,kernel_size=avr_width2) * sbig_spcal_day_coef
    img=ss.medfilt((img-masterdark.astype('float'))/masterflat.astype('float') - dark,kernel_size=avr_width2) * 0.045

    hdu_light = fits.PrimaryHDU(img)
    hdu_light.header['DATE-OBS']=f_date_iso
    hdu_light.header['BITPIX']=-64
    hdu_light.header['EXPTIME']=f_exp
    hdu_light.header['PROJ-TYP']='TAN'
    hdu_light.header['AMCAL-A0']=az0
    hdu_light.header['AMCAL-H0']=alt0
    hdu_light.header['AMCAL-A']=str(a)
    hdu_light.header['AMCAL-B']=str(b)
    hdu_light.header['AMCAL-C']=str(c)
    hdu_light.header['AMCAL-D']=str(d)
    hdu_light.header['SPCAL-DC']=sbig_spcal_day_coef
    hdu_light.header['SPCAL-DS']=sbig_spcal_std
    hdu_light.header['SPCAL-C']=sp_coef
    hdu_light.header['SPCAL-O']=sp_offset
    hdu_light.header['MED-WID1']=avr_width1
    hdu_light.header['MED-WID2']=avr_width2
    for j in range(len(pumping_scheme_list)):
        hdu_light.header['P-SCH-'+str(j)]=pumping_scheme_list[j]

    hdulist_light = fits.HDUList([hdu_light])
    hdulist_light.writeto(spath+fn.split('/')[-1].split('.')[-2]+".fit",overwrite=True)

    # plotting
    fig=plt.figure(figsize=(12.8,7.2))
    fig.set_size_inches(12.8, 7.2)
    dx=0.05
    ax1=plt.axes(position=[0.035000000000000003+dx/2, 0.26, 0.4, 0.7])
    pcm1=plt.pcolormesh(img,vmin=-20,vmax=20)
    ax1.set_ylim(img.shape[0],0)
    ax1.set_xlim(0,img.shape[1])
    plt.axis('equal')
    plt.title('LIGHT',loc='left')
    plt.title(f_date_iso,loc='right')

    ax1_cb=plt.axes(position=[0.035000000000000003+dx/2+0.4+0.008, 0.26, 0.015, 0.7])
    plt.colorbar(pcm1,ax1_cb)

    ax2=plt.axes(position=[0.485+0.035000000000000003+dx/2, 0.26, 0.4, 0.7])
    med=np.median(dark*sbig_spcal_day_coef)
    pcm2=plt.pcolormesh(dark*sbig_spcal_day_coef,vmin=med-100,vmax=med+100)
    
    ax2.set_ylim(img.shape[0],0)
    ax2.set_xlim(0,img.shape[1])
    plt.axis('equal')
    coefs_title=' '.join(("{:4.2f}".format(sbig_spcal_day_coef),
                          "{:5.3f}".format(sbig_spcal_std),
                          "{:4.2f}".format(sp_coef),
                          "{:4.2f}".format(sp_offset)))
    plt.title("SUBSTRACT",loc='left')
    plt.title("SP_CAL: " + coefs_title + "; MED=" + str(int(med)),loc='right')

    ax2_cb=plt.axes(position=[0.485+0.035000000000000003+dx/2+0.4+0.008, 0.26, 0.015, 0.7])
    plt.colorbar(pcm2,ax2_cb)

    ax3=plt.axes(position=[0.035000000000000003+dx/2, 0.05, 0.907, 0.17])

    xlim_left=f_x_date_start-12/60/24
    xlim_right=f_x_date_start+12/60/24
    ax3.set_xlim(xlim_left,xlim_right)
    ax3.set_ylim(0,1.)
    ax3.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
    ax3.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
    xfmt = dates.DateFormatter('%H:%M')
    ax3.xaxis.set_major_formatter(xfmt)
    for i in range(len(clean_pumping)):
        plt.plot(x_dates_clean[i],clean_pumping[i]*0.9,'r')

    plt.plot([f_x_date_start, f_x_date_end],[0.95,0.95],'b',lw=3)
    plt.grid()
    plt.gca().yaxis.set_major_locator(plt.NullLocator())

    plt.savefig(spath+fn.split('/')[-1].split('.')[-2]+".png")
    plt.close()

