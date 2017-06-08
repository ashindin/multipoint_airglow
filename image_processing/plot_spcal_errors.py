
# coding: utf-8

# In[1]:

import datetime
import numpy as np
from matplotlib import dates
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits


# In[2]:

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)


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
    fid_fit=fits.open(filename,ignore_missing_end=True);
    date_obs_str=fid_fit[0].header["DATE-OBS"]
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S")

def sbig_get_date_obs(filename):
    fid_fit=fits.open(filename);
    date_obs_str=fid_fit[0].header["DATE-OBS"]
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S.%f")

def s1c_get_exp_sec(filename):
    fdt=open(filename,'r',errors = 'ignore')
    fdt.seek(7*80)
    temp=fdt.read(80)
    exp_sec=float(temp.split("=")[-1].split(" ")[1])
    fdt.close()
    return exp_sec

def keo_get_exp_sec(filename):
    fid_fit=fits.open(filename,ignore_missing_end=True);
    exp_sec=float(fid_fit[0].header["EXPTIME"])
    fid_fit.close()
    return exp_sec

def sbig_get_exp_sec(filename):
    fid_fit=fits.open(filename);
    exp_sec=float(fid_fit[0].header["EXPTIME"])
    fid_fit.close()
    return exp_sec

def get_spcal_day_coefs(spcal_day_fname):
    fid = open(spcal_day_fname,'r')
    lines=fid.readlines()
    fid.close()
    lines_list=lines[1].split(" ")
    spcal_coef=float(lines_list[0])
    spcal_std=float(lines_list[1])
    return spcal_coef, spcal_std

def get_spcal_coefs(spcal_day_fname):
    fid = open(spcal_day_fname,'r')
    lines=fid.readlines()
    fid.close()
    spcal_coefs=[]
    for i in range(1,len(lines)):
        lines_list=lines[i].split(" ")
        spcal_coefs.append((lines_list[0],float(lines_list[1])))
    return spcal_coefs


# In[4]:

s1c_fit_path="../data/140824/s1c"
s1c_spcal_fname = "../spectrophotometric_calibration/s1c_140824.spcal"
s1c_spcal_day_fname = "../spectrophotometric_calibration/s1c_140824_day.spcal"
s1c_spcal_day_coef, s1c_spcal_std=get_spcal_day_coefs(s1c_spcal_day_fname)    
s1c_spcal_coefs=get_spcal_coefs(s1c_spcal_fname)
s1c_fit_filenames=[s1c_fit_path+"/"+s1c_spcal_coefs[i][0] for i in range(len(s1c_spcal_coefs))]
s1c_exp_sec=[s1c_get_exp_sec(f) for f in s1c_fit_filenames]
s1c_date_obs=[s1c_get_date_obs(s1c_fit_filenames[i],ut_shift=-4)+datetime.timedelta(seconds=s1c_exp_sec[i]/2) for i in range(len(s1c_spcal_coefs))]
s1c_x_dates=dates.date2num(s1c_date_obs)
s1c_sp_coefs=[sp[1] for sp in s1c_spcal_coefs]


# In[5]:

fig=plt.figure(figsize=(12.8,7.2))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
plt.plot(s1c_x_dates,s1c_sp_coefs,'ro')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef,s1c_spcal_day_coef],'k',lw=3,label='mean')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef-s1c_spcal_std,s1c_spcal_day_coef-s1c_spcal_std],'b',lw=3,ls="--",label='std')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef+s1c_spcal_std,s1c_spcal_day_coef+s1c_spcal_std],'b',lw=3,ls="--")

plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef-2*s1c_spcal_std,s1c_spcal_day_coef-2*s1c_spcal_std],'g',lw=3,ls="--",label='2std')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef+2*s1c_spcal_std,s1c_spcal_day_coef+2*s1c_spcal_std],'g',lw=3,ls="--")

plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef-3*s1c_spcal_std,s1c_spcal_day_coef-3*s1c_spcal_std],'gray',lw=3,ls="--",label='3std')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef+3*s1c_spcal_std,s1c_spcal_day_coef+3*s1c_spcal_std],'gray',lw=3,ls="--")


plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2014,8,24,17,43))
xlim_right=dates.date2num(datetime.datetime(2014,8,24,19,52))
ax.set_xlim(xlim_left,xlim_right)
plt.title("S1C 14/08/24 spectrophotometric calibration coefficients errors")
plt.ylabel("Calibration coefficient [R/ADCu]")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.grid()
# plt.show()
plt.savefig("s1c_140824_spcal_errors.png")
plt.close()


# In[6]:

keo_fit_path="../data/140824/keo"
keo_spcal_fname = "../spectrophotometric_calibration/keo_140824.spcal"
keo_spcal_day_fname = "../spectrophotometric_calibration/keo_140824_day.spcal"
keo_spcal_day_coef, keo_spcal_std=get_spcal_day_coefs(keo_spcal_day_fname)    
keo_spcal_coefs=get_spcal_coefs(keo_spcal_fname)
keo_fit_filenames=[keo_fit_path+"/"+keo_spcal_coefs[i][0] for i in range(len(keo_spcal_coefs))]
keo_exp_sec=[keo_get_exp_sec(f) for f in keo_fit_filenames]
keo_date_obs=[keo_get_date_obs(keo_fit_filenames[i])+datetime.timedelta(seconds=keo_exp_sec[i]/2) for i in range(len(keo_spcal_coefs))]
keo_x_dates=dates.date2num(keo_date_obs)
keo_sp_coefs=[sp[1] for sp in keo_spcal_coefs]


# In[7]:

fig=plt.figure(figsize=(12.8,7.2))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
plt.plot(keo_x_dates,keo_sp_coefs,'ro')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef,keo_spcal_day_coef],'k',lw=3,label='mean')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef-keo_spcal_std,keo_spcal_day_coef-keo_spcal_std],'b',lw=3,ls="--",label='std')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef+keo_spcal_std,keo_spcal_day_coef+keo_spcal_std],'b',lw=3,ls="--")

plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef-2*keo_spcal_std,keo_spcal_day_coef-2*keo_spcal_std],'g',lw=3,ls="--",label='2std')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef+2*keo_spcal_std,keo_spcal_day_coef+2*keo_spcal_std],'g',lw=3,ls="--")

plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef-3*keo_spcal_std,keo_spcal_day_coef-3*keo_spcal_std],'gray',lw=3,ls="--",label='3std')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef+3*keo_spcal_std,keo_spcal_day_coef+3*keo_spcal_std],'gray',lw=3,ls="--")


plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2014,8,24,17,16))
xlim_right=dates.date2num(datetime.datetime(2014,8,24,20,43))
ax.set_xlim(xlim_left,xlim_right)
plt.title("KEO 14/08/24 spectrophotometric calibration coefficients errors")
plt.ylabel("Calibration coefficient [R/ADCu]")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.grid()
# plt.show()
plt.savefig("keo_140824_spcal_errors.png")
plt.close()


# In[8]:

s1c_fit_path="../data/140826/s1c"
s1c_spcal_fname = "../spectrophotometric_calibration/s1c_140826.spcal"
s1c_spcal_day_fname = "../spectrophotometric_calibration/s1c_140826_day.spcal"
s1c_spcal_day_coef, s1c_spcal_std=get_spcal_day_coefs(s1c_spcal_day_fname)    
s1c_spcal_coefs=get_spcal_coefs(s1c_spcal_fname)
s1c_fit_filenames=[s1c_fit_path+"/"+s1c_spcal_coefs[i][0] for i in range(len(s1c_spcal_coefs))]
s1c_exp_sec=[s1c_get_exp_sec(f) for f in s1c_fit_filenames]
s1c_date_obs=[s1c_get_date_obs(s1c_fit_filenames[i],ut_shift=-4)+datetime.timedelta(seconds=s1c_exp_sec[i]/2) for i in range(len(s1c_spcal_coefs))]
s1c_x_dates=dates.date2num(s1c_date_obs)
s1c_sp_coefs=[sp[1] for sp in s1c_spcal_coefs]


# In[9]:

fig=plt.figure(figsize=(12.8,7.2))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
plt.plot(s1c_x_dates,s1c_sp_coefs,'ro')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef,s1c_spcal_day_coef],'k',lw=3,label='mean')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef-s1c_spcal_std,s1c_spcal_day_coef-s1c_spcal_std],'b',lw=3,ls="--",label='std')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef+s1c_spcal_std,s1c_spcal_day_coef+s1c_spcal_std],'b',lw=3,ls="--")

plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef-2*s1c_spcal_std,s1c_spcal_day_coef-2*s1c_spcal_std],'g',lw=3,ls="--",label='2std')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef+2*s1c_spcal_std,s1c_spcal_day_coef+2*s1c_spcal_std],'g',lw=3,ls="--")

plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef-3*s1c_spcal_std,s1c_spcal_day_coef-3*s1c_spcal_std],'gray',lw=3,ls="--",label='3std')
plt.plot([s1c_x_dates[0],s1c_x_dates[-1]],[s1c_spcal_day_coef+3*s1c_spcal_std,s1c_spcal_day_coef+3*s1c_spcal_std],'gray',lw=3,ls="--")


plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2014,8,26,18,10))
xlim_right=dates.date2num(datetime.datetime(2014,8,26,21,19))
ax.set_xlim(xlim_left,xlim_right)
ax.set_ylim(2.0,2.7)
plt.title("S1C 14/08/26 spectrophotometric calibration coefficients errors")
plt.ylabel("Calibration coefficient [R/ADCu]")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.grid()
# plt.show()
plt.savefig("s1c_140826_spcal_errors.png")
plt.close()


# In[10]:

keo_fit_path="../data/140826/keo"
keo_spcal_fname = "../spectrophotometric_calibration/keo_140826.spcal"
keo_spcal_day_fname = "../spectrophotometric_calibration/keo_140826_day.spcal"
keo_spcal_day_coef, keo_spcal_std=get_spcal_day_coefs(keo_spcal_day_fname)    
keo_spcal_coefs=get_spcal_coefs(keo_spcal_fname)
keo_fit_filenames=[keo_fit_path+"/"+keo_spcal_coefs[i][0] for i in range(len(keo_spcal_coefs))]
keo_exp_sec=[keo_get_exp_sec(f) for f in keo_fit_filenames]
keo_date_obs=[keo_get_date_obs(keo_fit_filenames[i])+datetime.timedelta(seconds=keo_exp_sec[i]/2) for i in range(len(keo_spcal_coefs))]
keo_x_dates=dates.date2num(keo_date_obs)
keo_sp_coefs=[sp[1] for sp in keo_spcal_coefs]


# In[11]:

fig=plt.figure(figsize=(12.8,7.2))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
plt.plot(keo_x_dates,keo_sp_coefs,'ro')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef,keo_spcal_day_coef],'k',lw=3,label='mean')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef-keo_spcal_std,keo_spcal_day_coef-keo_spcal_std],'b',lw=3,ls="--",label='std')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef+keo_spcal_std,keo_spcal_day_coef+keo_spcal_std],'b',lw=3,ls="--")

plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef-2*keo_spcal_std,keo_spcal_day_coef-2*keo_spcal_std],'g',lw=3,ls="--",label='2std')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef+2*keo_spcal_std,keo_spcal_day_coef+2*keo_spcal_std],'g',lw=3,ls="--")

plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef-3*keo_spcal_std,keo_spcal_day_coef-3*keo_spcal_std],'gray',lw=3,ls="--",label='3std')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef+3*keo_spcal_std,keo_spcal_day_coef+3*keo_spcal_std],'gray',lw=3,ls="--")


plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2014,8,26,17,29))
xlim_right=dates.date2num(datetime.datetime(2014,8,26,21,49))
ax.set_xlim(xlim_left,xlim_right)
ax.set_ylim(0.3,0.7)
plt.title("KEO 14/08/26 spectrophotometric calibration coefficients errors")
plt.ylabel("Calibration coefficient [R/ADCu]")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.grid()
# plt.show()
plt.savefig("keo_140826_spcal_errors.png")
plt.close()


# In[12]:

keo_fit_path="../data/160829/keo"
keo_spcal_fname = "../spectrophotometric_calibration/keo_160829.spcal"
keo_spcal_day_fname = "../spectrophotometric_calibration/keo_160829_day.spcal"
keo_spcal_day_coef, keo_spcal_std=get_spcal_day_coefs(keo_spcal_day_fname)    
keo_spcal_coefs=get_spcal_coefs(keo_spcal_fname)
keo_fit_filenames=[keo_fit_path+"/"+keo_spcal_coefs[i][0] for i in range(len(keo_spcal_coefs))]
keo_exp_sec=[keo_get_exp_sec(f) for f in keo_fit_filenames]
keo_date_obs=[keo_get_date_obs(keo_fit_filenames[i])+datetime.timedelta(seconds=keo_exp_sec[i]/2) for i in range(len(keo_spcal_coefs))]
keo_x_dates=dates.date2num(keo_date_obs)
keo_sp_coefs=[sp[1] for sp in keo_spcal_coefs]


# In[13]:

fig=plt.figure(figsize=(12.8,7.2))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
plt.plot(keo_x_dates,keo_sp_coefs,'ro')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef,keo_spcal_day_coef],'k',lw=3,label='mean')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef-keo_spcal_std,keo_spcal_day_coef-keo_spcal_std],'b',lw=3,ls="--",label='std')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef+keo_spcal_std,keo_spcal_day_coef+keo_spcal_std],'b',lw=3,ls="--")

plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef-2*keo_spcal_std,keo_spcal_day_coef-2*keo_spcal_std],'g',lw=3,ls="--",label='2std')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef+2*keo_spcal_std,keo_spcal_day_coef+2*keo_spcal_std],'g',lw=3,ls="--")

plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef-3*keo_spcal_std,keo_spcal_day_coef-3*keo_spcal_std],'gray',lw=3,ls="--",label='3std')
plt.plot([keo_x_dates[0],keo_x_dates[-1]],[keo_spcal_day_coef+3*keo_spcal_std,keo_spcal_day_coef+3*keo_spcal_std],'gray',lw=3,ls="--")


plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2016,8,29,16,40))
xlim_right=dates.date2num(datetime.datetime(2016,8,29,21,28))
ax.set_xlim(xlim_left,xlim_right)
ax.set_ylim(0.4,0.8)
plt.title("KEO 16/08/29 spectrophotometric calibration coefficients errors")
plt.ylabel("Calibration coefficient [R/ADCu]")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.grid()
# plt.show()
plt.savefig("keo_160829_spcal_errors.png")
plt.close()


# In[14]:

sbig_fit_path="../data/160829/sbig"
sbig_spcal_fname = "../spectrophotometric_calibration/sbig_160829.spcal"
sbig_spcal_day_fname = "../spectrophotometric_calibration/sbig_160829_day.spcal"
sbig_spcal_day_coef, sbig_spcal_std=get_spcal_day_coefs(sbig_spcal_day_fname)    
sbig_spcal_coefs=get_spcal_coefs(sbig_spcal_fname)
sbig_fit_filenames=[sbig_fit_path+"/"+sbig_spcal_coefs[i][0] for i in range(len(sbig_spcal_coefs))]
sbig_exp_sec=[sbig_get_exp_sec(f) for f in sbig_fit_filenames]
sbig_date_obs=[sbig_get_date_obs(sbig_fit_filenames[i])+datetime.timedelta(seconds=sbig_exp_sec[i]/2) for i in range(len(sbig_spcal_coefs))]
sbig_x_dates=dates.date2num(sbig_date_obs)
sbig_sp_coefs=[sp[1] for sp in sbig_spcal_coefs]


# In[15]:

fig=plt.figure(figsize=(12.8,7.2))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
plt.plot(sbig_x_dates,sbig_sp_coefs,'ro')
plt.plot([sbig_x_dates[0],sbig_x_dates[-1]],[sbig_spcal_day_coef,sbig_spcal_day_coef],'k',lw=3,label='mean')
plt.plot([sbig_x_dates[0],sbig_x_dates[-1]],[sbig_spcal_day_coef-sbig_spcal_std,sbig_spcal_day_coef-sbig_spcal_std],'b',lw=3,ls="--",label='std')
plt.plot([sbig_x_dates[0],sbig_x_dates[-1]],[sbig_spcal_day_coef+sbig_spcal_std,sbig_spcal_day_coef+sbig_spcal_std],'b',lw=3,ls="--")

plt.plot([sbig_x_dates[0],sbig_x_dates[-1]],[sbig_spcal_day_coef-2*sbig_spcal_std,sbig_spcal_day_coef-2*sbig_spcal_std],'g',lw=3,ls="--",label='2std')
plt.plot([sbig_x_dates[0],sbig_x_dates[-1]],[sbig_spcal_day_coef+2*sbig_spcal_std,sbig_spcal_day_coef+2*sbig_spcal_std],'g',lw=3,ls="--")

plt.plot([sbig_x_dates[0],sbig_x_dates[-1]],[sbig_spcal_day_coef-3*sbig_spcal_std,sbig_spcal_day_coef-3*sbig_spcal_std],'gray',lw=3,ls="--",label='3std')
plt.plot([sbig_x_dates[0],sbig_x_dates[-1]],[sbig_spcal_day_coef+3*sbig_spcal_std,sbig_spcal_day_coef+3*sbig_spcal_std],'gray',lw=3,ls="--")


plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2016,8,29,17,13))
xlim_right=dates.date2num(datetime.datetime(2016,8,29,20,22))
ax.set_xlim(xlim_left,xlim_right)
ax.set_ylim(0.035,0.055)
plt.title("SBIG 16/08/29 spectrophotometric calibration coefficients errors")
plt.ylabel("Calibration coefficient [R/ADCu]")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.grid()
# plt.show()
plt.savefig("sbig_160829_spcal_errors.png")
plt.close()

