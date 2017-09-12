
# coding: utf-8

# In[2]:

import datetime
import os, sys, inspect
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.signal as ss
from matplotlib import dates

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../astrometric_calibration")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from tan_module import *
# from sp_cal_module import *


# In[15]:

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


# In[3]:

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
def keo_get_date_obs(filename):
    fid_fit=fits.open(filename);
    date_obs_str=fid_fit[0].header["DATE-OBS"]
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S")
def sbig_get_date_obs(filename):
    fid_fit=fits.open(filename);
    date_obs_str=fid_fit[0].header["DATE-OBS"]
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S.%f")
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


# In[18]:

pumping_scheme_file="../image_processing/pump140824.scheme"
pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)
date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)


# In[78]:

def figure4(fit_path,masterdark_fname,masterflat_fname,R,xpix=144,ypix=144,rang=15):
    fit_filenames=[fit_path+'/'+fn for fn in next(os.walk(fit_path))[2]]
#     fid=open(solve_pars_fname,'r')
#     lines=fid.readlines()[1::]
#     fid.close()
#     lines=[line[:-1] for line in lines]
    masterdark=fits.getdata(masterdark_fname,ignore_missing_end=True).astype('float')
    masterflat=fits.getdata(masterflat_fname,ignore_missing_end=True).astype('float')
    
#     cp_ind=[35+i*24 for i in range(16)]+[443,514]
    cp_ind=[35, 58, 83, 107, 130, 155, 175, 203, 227, 251, 275, 299, 323, 347, 371,
            390, 443, 498]
    
    Int=[]
    dInt=[]
    xdates=[]
    cp_x=[]
    cp_int=[]
    back_int=[]

    for i in range(len(fit_filenames)):
        fit_fn=fit_filenames[i]
        dateobs=s1c_get_date_obs(fit_fn,ut_shift=-4)
        img=fits.getdata(fit_fn,ignore_missing_end=True).astype(float)
        img=(img-masterdark)/masterflat*R
        xdate=dates.date2num(dateobs+datetime.timedelta(seconds=7.5))
        xdates.append(xdate)
        I=np.median(img[ypix-rang:ypix+rang+1,xpix+rang+1])
        Int.append(I)
        if i in cp_ind:
            cp_x.append(xdate)
            cp_int.append(I)     
    
    for i in range(len(Int)):
        back_int.append(np.interp(xdates[i],cp_x,cp_int))
        dInt.append(Int[i]-back_int[i])
    
    hei=0.36
    fig=plt.figure(figsize=(6.69,6))
    ax=plt.axes(position=[0.092, 0.05+0.06+hei+0.125, 0.88, hei])    
    ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(1,61,6)))
    ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(4,64,6)))
    xfmt = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)    
    plt.plot(xdates,Int,'ro',ms=2)
    plt.plot(cp_x,cp_int,'b')
    plt.plot(cp_x,cp_int,'bo',ms=4)
#     plt.plot(xdates,back_int,'b')
    xlim_left=dates.date2num(datetime.datetime(2014,8,24,17,50))
    xlim_right=dates.date2num(datetime.datetime(2014,8,24,19,46))    
    plt.ylim([400,650])
    ax.set_xlim(xlim_left,xlim_right)
    plt.ylabel("$B$, Рэлей")
#     plt.xlabel("Время ЧЧ:ММ, UTC")
    plt.xticks(rotation=90)
    for i in range(len(clean_pumping)):
        plt.plot(x_dates_clean[i],clean_pumping[i]*25+400,'k')
    plt.grid()
    ax.annotate('а',xy=(.05, .98), xycoords='figure fraction')

    ax2=plt.axes(position=[0.092, 0.05+0.08, 0.88, hei])    
    ax2.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(1,61,6)))
    ax2.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(4,64,6)))
    xfmt = dates.DateFormatter('%H:%M')
    ax2.xaxis.set_major_formatter(xfmt)    
    plt.plot(xdates,dInt,'r')
#     plt.plot(cp_x,cp_int,'b')
    xlim_left=dates.date2num(datetime.datetime(2014,8,24,17,50))
    xlim_right=dates.date2num(datetime.datetime(2014,8,24,19,46))    
    plt.ylim([-25,45])
    ax2.set_xlim(xlim_left,xlim_right)
    plt.ylabel("$\Delta B$, Рэлей")
    plt.xlabel("Время ЧЧ:ММ, UTC")
    plt.xticks(rotation=90)
    for i in range(len(clean_pumping)):
        plt.plot(x_dates_clean[i],clean_pumping[i]*10-25,'k')
    plt.grid()
    ax2.annotate('б',xy=(.05, .50), xycoords='figure fraction')
    plt.savefig("fig4.png",dpi=300)
    plt.savefig("fig4.eps",dpi=300)
    plt.savefig("fig4.pdf",dpi=300)
#     plt.show()
#     print(cp_x)
#     print(cp_int)
    return cp_x, cp_int
fit_path="../data/140824/s1c"
masterdark_fname="../spectrophotometric_calibration/s1c_140824_masterdark.fit"
masterflat_fname="../spectrophotometric_calibration/s1c_master.flat"
R=2.20
# solve_pars_fname="../astrometric_calibration/s1c_140824_solve.pars"
# save_fname="../spectrophotometric_calibration/s1c_140824.spcal"
cp_x, cp_int = figure4(fit_path,masterdark_fname,masterflat_fname,R)
# 35,59,83
