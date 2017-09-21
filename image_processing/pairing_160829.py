
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
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import dates


# In[2]:

def get_fit_pars(fname):
#     print(fname)
    fid_fit=fits.open(fname);
    date_obs=datetime.datetime.strptime(fid_fit[0].header["DATE-OBS"],"%Y-%m-%dT%H:%M:%S.%f")
#     print(date_obs)
    exp_sec=fid_fit[0].header["EXPTIME"]
    sp_offset=fid_fit[0].header["SPCAL-O"]
    fid_fit.close()    
    return date_obs, exp_sec, sp_offset
    
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

def make_couple_fit(couple_fn, fn1, fn21):
    if os.path.isfile(couple_fn):
        os.remove(couple_fn)
    
    fn1d=fits.getdata(fn1)
    fn1h=fits.getheader(fn1)
    fn1h['LAT-DEG']=55.9305361
    fn1h['LON-DEG']=48.7444861
    fn1h['HEI_M']=91.
    fn1_date=datetime.datetime.strptime(fn1h['DATE-OBS'],'%Y-%m-%dT%H:%M:%S.%f')
    fn1h['DATE-OBS']=fn1_date.strftime('%Y-%m-%dT%H:%M:%S.%f')

    fn21d=fits.getdata(fn21)
    fn21h=fits.getheader(fn21)
    fn21_exptime=fn21h['EXPTIME']
    fn21_spcalc=fn21h['SPCAL-C']
    fn21_spcalo=fn21h['SPCAL-O']
    fn21_date_start=datetime.datetime.strptime(fn21h['DATE-OBS'],'%Y-%m-%dT%H:%M:%S.%f')
    fn21_date=fn21_date_start+datetime.timedelta(seconds=fn21_exptime/2)
    fn21_xdate=dates.date2num(fn21_date)
    fn21_date_end=fn21_date_start+datetime.timedelta(seconds=fn21_exptime)
    
    fn21h['LAT-DEG']=56.1501667
    fn21h['LON-DEG']=46.1050833
    fn21h['HEI_M']=183.
    
    fits.append(couple_fn,fn1d,fn1h)
    fits.append(couple_fn,fn21d,fn21h)  

    return None


# In[4]:

keo_fit_path="./keo160829_glowfit"
keo_fit_filenames=[]
keo_fit_fullnames=[]
for fn in next(os.walk(keo_fit_path))[2]:
    if fn[-3:]=='fit' and fn[-9:]!='_dark.fit':
        keo_fit_filenames.append(fn)
        keo_fit_fullnames.append(keo_fit_path + '/' + fn)
keo_fit_filenames=sorted(keo_fit_filenames)
keo_fit_fullnames=sorted(keo_fit_fullnames)

# sorted([keo_fit_path+'/'+fn for fn in next(os.walk(keo_fit_path))[2]])
sbig_fit_path="./sbig160829_glowfit"
sbig_fit_filenames=[]
sbig_fit_fullnames=[]
for fn in next(os.walk(sbig_fit_path))[2]:
    if fn[-3:]=='fit' and fn[-9:]!='_dark.fit':
        sbig_fit_filenames.append(fn)
        sbig_fit_fullnames.append(sbig_fit_path + '/' + fn)
sbig_fit_filenames=sorted(sbig_fit_filenames)
sbig_fit_fullnames=sorted(sbig_fit_fullnames)

# In[5]:

spath="./couples160829/"
if not os.path.exists(spath):
    os.makedirs(spath)


# In[6]:

sbig_date_start=[]
sbig_date=[]
sbig_x_date=[]
sbig_date_end=[]
sbig_sp_offsets=[]
for i in range(len(sbig_fit_filenames)):
    fn = sbig_fit_filenames[i]
    date_obs_str, exp_sec, sp_offset = get_fit_pars(sbig_fit_path + '/' + fn)
#     print(date_obs_str)
    sbig_date_start.append(date_obs_str)
    sbig_date.append(date_obs_str+datetime.timedelta(seconds=exp_sec/2))
    sbig_x_date.append(dates.date2num(sbig_date[-1]))
    sbig_date_end.append(date_obs_str+datetime.timedelta(seconds=exp_sec))
    sbig_sp_offsets.append(sp_offset)
sbig_x_date=np.array(sbig_x_date)
#     print(sbig_sp_offsets[-1],sbig_date_start[-1],sbig_date[-1],sbig_date_end[-1])


# In[ ]:

c_i=0
for i in range(len(keo_fit_filenames)):
# for i in range(60,61):
    fn = keo_fit_filenames[i]
#     print(fn)
    date_obs_str, exp_sec, sp_offset = get_fit_pars(keo_fit_path + '/' + fn)
    if sp_offset>3:
        continue    
    keo_date=date_obs_str+datetime.timedelta(seconds=exp_sec/2)
    keo_x_date=dates.date2num(keo_date)
    sbig_dates_offset=np.abs(sbig_x_date-keo_x_date)
    sbig_arg_sort=np.argsort(sbig_dates_offset)
#     print(sbig_arg_sort)
    if sbig_sp_offsets[sbig_arg_sort[0]]>3 and sbig_sp_offsets[sbig_arg_sort[0]]<10:
        continue
    if sbig_dates_offset[sbig_arg_sort[0]]>10/60/60/24:
        continue
#     print(fn, sbig_fit_filenames[sbig_arg_sort[0]], sbig_fit_filenames[sbig_arg_sort[1]])
    c_i+=1
    couple_name=spath + 'couple_160829_' + '{:03}'.format(c_i) + '.fits'
#     print(couple_name)
#     print(keo_fit_fullnames[i],sbig_fit_fullnames[sbig_arg_sort[0]], sbig_fit_fullnames[sbig_arg_sort[1]])
    make_couple_fit(couple_name,keo_fit_fullnames[i],sbig_fit_fullnames[sbig_arg_sort[0]])
#     print(c_i)


# In[8]:

couples_fit_filenames=[]
for fn in next(os.walk(spath))[2]:
    if fn[-4:]=='fits':
        couples_fit_filenames.append(spath+fn)
couples_fit_filenames=sorted(couples_fit_filenames)


# In[ ]:

# plotting
vmin=-15
vmax=15
for i in range(len(couples_fit_filenames)):
    # for i in range(0,1):
    sys.stdout.write('\r')
    sys.stdout.write("Processing frame "+str(i+1)+"/"+str(len(couples_fit_filenames)))
    sys.stdout.flush()

    fn=couples_fit_filenames[i]
    img1=fits.getdata(fn,0)
    img2=fits.getdata(fn,1)
    
    dateobs1_str=fits.getval(fn,'DATE-OBS',0)
    dateobs2_str=fits.getval(fn,'DATE-OBS',1)
    dateobs1=datetime.datetime.strptime(dateobs1_str,"%Y-%m-%dT%H:%M:%S.%f")
    dateobs2=datetime.datetime.strptime(dateobs2_str,"%Y-%m-%dT%H:%M:%S.%f")
    dateobs1x=dates.date2num(dateobs1)
    dateobs2x=dates.date2num(dateobs2)
    
    spcalo1=fits.getval(fn,'SPCAL-O',0)
    spcalo2=fits.getval(fn,'SPCAL-O',1)
    
    exptime1=fits.getval(fn,'EXPTIME',0)
    exptime2=fits.getval(fn,'EXPTIME',1)
    dateend1=dateobs1+datetime.timedelta(seconds=exptime1)
    dateend2=dateobs2+datetime.timedelta(seconds=exptime2)
    dateend1x=dates.date2num(dateend1)
    dateend2x=dates.date2num(dateend2)
        
    pumping_scheme_str_list=[]
    s_i=0
    while True:
        try:
            temp=fits.getval(fn,'P-SCH-'+str(s_i),0)
        except:
            break
        pumping_scheme_str_list.append(temp)
        s_i+=1
    
    pumping_scheme=[]
    for i in range(len(pumping_scheme_str_list)):
        ps_str=pumping_scheme_str_list[i]
        ps_str_list=ps_str.split(',')
        ds=datetime.datetime.strptime(ps_str_list[0],"%Y-%m-%dT%H:%M:%S")
        tau=float(ps_str_list[1].split(' min')[0])
        T=float(ps_str_list[2].split(' min')[0])
        N=float(ps_str_list[3])
        V=float(ps_str_list[4])
        F=float(ps_str_list[5].split(' kHz')[0])
        pumping_scheme.append([ds,tau,T,N,V,F])
    
    date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
    
    fig=plt.figure(figsize=(12.8,7.2))
    fig.set_size_inches(12.8, 7.2)
    dx=0.05
    ax1=plt.axes(position=[0.035000000000000003+dx/2, 0.26, 0.4, 0.7])    
    pcm1=plt.pcolormesh(img1,vmin=vmin,vmax=vmax)
    plt.title("KEO: "+ dateobs1_str, loc='left')
    plt.title("SP_OFFSET: " + "{:4.2f}".format(spcalo1),loc='right')
#     ax1.set_ylim(img1.shape[0],0)
#     ax1.set_xlim(img1.shape[1],0)
    ax1.set_ylim(425,325)
    ax1.set_xlim(275,175)
#     plt.axis('equal')
    ax1_cb=plt.axes(position=[0.035000000000000003+dx/2+0.4+0.008, 0.26, 0.015, 0.7])
    plt.colorbar(pcm1,ax1_cb)
    
    ax2=plt.axes(position=[0.485+0.035000000000000003+dx/2, 0.26, 0.4, 0.7])
    pcm2=plt.pcolormesh(img2,vmin=vmin,vmax=vmax)
    plt.title("sbig: "+ dateobs2_str,loc='left')
    plt.title("SP_OFFSET: " + "{:4.2f}".format(spcalo2),loc='right')    
    ax2.set_ylim(img2.shape[0],0)
    ax2.set_xlim(0,img2.shape[1])
    plt.axis('equal')
    ax2_cb=plt.axes(position=[0.485+0.035000000000000003+dx/2+0.4+0.008, 0.26, 0.015, 0.7])
    plt.colorbar(pcm2,ax2_cb)    
    
    ax3=plt.axes(position=[0.035000000000000003+dx/2, 0.05, 0.907, 0.17])
    xlim_left=dateobs1x-12/60/24
    xlim_right=dateobs1x+12/60/24
    ax3.set_xlim(xlim_left,xlim_right)
    ax3.set_ylim(0,1.)
    ax3.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
    ax3.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
    xfmt = dates.DateFormatter('%H:%M')
    ax3.xaxis.set_major_formatter(xfmt)
    for i in range(len(clean_pumping)):
        plt.plot(x_dates_clean[i],clean_pumping[i]*0.9,'r')

    plt.plot([dateobs1x, dateend1x],[0.95,0.95],'b',lw=3)
    plt.plot([dateobs2x, dateend2x],[0.92,0.92],'g',lw=3)
    plt.grid()
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    
    plt.savefig(fn[0:-4]+"png")
    plt.close()    


# In[ ]:



