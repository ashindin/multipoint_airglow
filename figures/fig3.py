
# coding: utf-8

# In[1]:

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


# In[2]:

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


# In[3]:

def figure3(fit_path,solve_pars_fname):
    fid=open(solve_pars_fname,'r')
    lines=fid.readlines()[1::]
    fid.close()
    lines=[line[:-1] for line in lines]
    fig=plt.figure(figsize=(6.69,3.4))
    ax=plt.axes(position=[0.092, 0.1+0.12, 0.90, 0.75])    
    ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(2,62,6)))
#     ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
    xfmt = dates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    R=[]
    R_st=[]
    for line in lines:
        line_list=line.split()
        fit_fn=fit_path+"/"+line_list[0]
        R_coef=float(line_list[1])
        R.append(R_coef)
        R_std=float(line_list[2])
        R_st.append(R_std)
        dateobs=s1c_get_date_obs(fit_fn,ut_shift=-4)
        xdate=dates.date2num(dateobs)
        plt.plot(xdate,R_coef,'ro',ms=1)
        plt.plot([xdate,xdate],[R_coef-R_std, R_coef+R_std],'r',lw=0.5)
        dx=0.0003
        plt.plot([xdate-dx,xdate+dx],[R_coef-R_std, R_coef-R_std],'r',lw=0.5)
        plt.plot([xdate-dx,xdate+dx],[R_coef+R_std, R_coef+R_std],'r',lw=0.5)
    
    R_med=np.median(R[500::])
    R_st_med=np.median(R_st[500::])
    R_std2=np.std(R[500::])
    xlim_left=dates.date2num(datetime.datetime(2014,8,26,18,10))
    xlim_right=dates.date2num(datetime.datetime(2014,8,26,21,15))
    
    plt.plot([xlim_left,xlim_right],[R_med,R_med],'k',lw=1.5)
    plt.plot([xlim_left,xlim_right],[R_med-R_std2-R_st_med,R_med-R_std2-R_st_med],'k',lw=1,ls=":")
    plt.plot([xlim_left,xlim_right],[R_med+R_std2+R_st_med,R_med+R_std2+R_st_med],'k',lw=1,ls=":")
    plt.ylim([2,4])

    ax.set_xlim(xlim_left,xlim_right)
    plt.ylabel("$R$")
    plt.xlabel("Время ЧЧ:ММ, UTC")
    plt.xticks(rotation=90)
    plt.grid()
    plt.savefig("fig3.png",dpi=300)
    plt.savefig("fig3.eps",dpi=300)
    plt.savefig("fig3.pdf",dpi=300)
    #plt.show()
    return R_med, R_std2+R_st_med
fit_path="../data/140826/s1c"
solve_pars_fname="../astrometric_calibration/s1c_140826_solve.pars"
save_fname="../spectrophotometric_calibration/s1c_140826.spcal"
figure3(fit_path,save_fname)    


# In[4]:

# R_median=R_median[500::]
# R_filt=R_median[np.where(R_median>0)]
# R_day=np.median(R_filt)
# R_std2=np.std(R_filt)
# print(R_day, R_std2)
# fid=open(save_fname,'w')
# fid.write("# Median camera calibration coefficient [Rayleighs per ADC unit] and its std for 14/08/26:\n")
# fid.write(str(R_day)+" "+str(R_std2))
# fid.close()

