
# coding: utf-8

# In[1]:

import datetime, os, copy
import numpy as np


# In[2]:

from astropy.io import fits


# In[3]:

from matplotlib import dates
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# In[4]:

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

def get_mod_pumping_scheme(pumping_scheme,prohibited_area_min):
    mod_pumping_scheme=copy.deepcopy(pumping_scheme)
    for i in range(len(mod_pumping_scheme)):
        mod_pumping_scheme[i][1]+=prohibited_area_min
    return mod_pumping_scheme

def get_frame_time_axe(date_obs,exp_sec):
    if date_obs.microsecond>500000:
        st_datetime=datetime.datetime(date_obs.year,date_obs.month,date_obs.day,date_obs.hour,date_obs.minute,date_obs.second+1)
    else:
        st_datetime=datetime.datetime(date_obs.year,date_obs.month,date_obs.day,date_obs.hour,date_obs.minute,date_obs.second)
    date_axe=[st_datetime+datetime.timedelta(seconds=i) for i in range(int(exp_sec))]
    x_dates=dates.date2num(date_axe)
    return date_axe, x_dates

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

def get_overlap_coef(frame_date_axe, date_axe, pumping):
    overlap_coef=0.
    for dt in frame_date_axe:
        for i in range(len(pumping)):
            for j in range(len(date_axe[i])):
                if dt==date_axe[i][j]:
                    overlap_coef+=pumping[i][j]
    overlap_coef/=len(frame_date_axe)
    return overlap_coef

def get_position_on_scheme(dt,pumping_scheme):
    ps_st_x_dates=np.zeros(len(pumping_scheme))
    dtx=dates.date2num(dt)
    ps_index=-1
    period_index=-1
    period_offset=0.
    
    for i in range(len(pumping_scheme)):
        ps_st_x_dates[i]=dates.date2num(pumping_scheme[i][0])
        if dtx>ps_st_x_dates[i]:
            ps_index=i
    if ps_index>=0:
        ps=pumping_scheme[ps_index]
        T=ps[2]
        num_p=int(ps[3])
        period_x_dates=np.zeros(num_p)
        for i in range(num_p):
            period_x_dates[i]=dates.date2num(ps[0]+datetime.timedelta(seconds=i*T*60))
            if dtx>period_x_dates[i]:
                period_index=i                
        if period_index>=0:
            period_offset=(dtx-period_x_dates[period_index])/(T*1./24/60)
    return ps_index, period_index, period_offset


# # 14/08/24

# In[5]:

# aug 24/08/14
pumping_scheme_file="pump140824.scheme"
s1c_fit_path="../data/140824/s1c"
keo_fit_path="../data/140824/keo"

s1c_spcal_fname = "../spectrophotometric_calibration/s1c_140824.spcal"
s1c_spcal_day_fname = "../spectrophotometric_calibration/s1c_140824_day.spcal"
s1c_spcal_day_coef, s1c_spcal_std=get_spcal_day_coefs(s1c_spcal_day_fname)    
s1c_spcal_coefs=get_spcal_coefs(s1c_spcal_fname)   

keo_spcal_fname = "../spectrophotometric_calibration/keo_140824.spcal"
keo_spcal_day_fname = "../spectrophotometric_calibration/keo_140824_day.spcal"
keo_spcal_day_coef, keo_spcal_std=get_spcal_day_coefs(keo_spcal_day_fname)
keo_spcal_coefs=get_spcal_coefs(keo_spcal_fname)

prohibited_area_min=2.
pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)
mod_pumping_scheme=get_mod_pumping_scheme(pumping_scheme,prohibited_area_min)
# pumping_scheme, mod_pumping_scheme


# In[6]:

s1c_fit_filenames=sorted([s1c_fit_path+'/'+fn for fn in next(os.walk(s1c_fit_path))[2]])
keo_fit_filenames=sorted([keo_fit_path+'/'+fn for fn in next(os.walk(keo_fit_path))[2]])

s1c_date_obs=[s1c_get_date_obs(f,ut_shift=-4) for f in s1c_fit_filenames]
s1c_exp_sec=[s1c_get_exp_sec(f) for f in s1c_fit_filenames]

keo_date_obs=[keo_get_date_obs(f) for f in keo_fit_filenames]
keo_exp_sec=[keo_get_exp_sec(f) for f in keo_fit_filenames]

date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
mod_date_axe_clean, mod_x_dates_clean, mod_clean_pumping = make_clean_pumping_scheme(mod_pumping_scheme)


# In[7]:

sum_of_periods=int(sum([ps[3] for ps in mod_pumping_scheme]))
BF_s1c=["" for i in range(sum_of_periods+1)]

for frame_ind in range(len(s1c_fit_filenames)):
    fn=s1c_fit_filenames[frame_ind].split("/")[-1]
    #
    
    ps_index, period_index, period_offset = get_position_on_scheme(s1c_date_obs[frame_ind]+datetime.timedelta(seconds=s1c_exp_sec[frame_ind]/2),mod_pumping_scheme)
    if ps_index==-1:
        BF_s1c_ind=0
    else:
        BF_s1c_ind=1
#         for i in range(ps_index):
        BF_s1c_ind+=int(sum([ps[3] for ps in mod_pumping_scheme[0:ps_index]]))
        BF_s1c_ind+=period_index
    #
    frame_date_axe, frame_x_dates = get_frame_time_axe(s1c_date_obs[frame_ind], s1c_exp_sec[frame_ind])
    overlap_coef=get_overlap_coef(frame_date_axe, mod_date_axe_clean, mod_clean_pumping)
    if overlap_coef<0.33 and period_offset<=1.:
        for i in range(len(s1c_spcal_coefs)):
            if s1c_spcal_coefs[i][0]==fn:
                i_find=i
                break
#                 if frame_ind==82: print("WoW!",abs(s1c_spcal_coefs[i][1]-s1c_spcal_day_coef),2*s1c_spcal_std);
        if abs(s1c_spcal_coefs[i_find][1]-s1c_spcal_day_coef)<=3*s1c_spcal_std:
            BF_s1c[BF_s1c_ind]=fn
            sign_flag=True
#     if i_find>-1: print(frame_ind,fn,s1c_date_obs[frame_ind],overlap_coef,period_offset,sign_flag,abs(s1c_spcal_coefs[i_find][1]-s1c_spcal_day_coef),2*s1c_spcal_std);
    i_find=-1
    sign_flag=False


# In[8]:

# BF_s1c
save_filename='s1c_140824_base.frames'
fid=open(save_filename,'w')
for bf in BF_s1c:
    fid.write(bf+'\n')
fid.close()


# In[9]:

sum_of_periods=int(sum([ps[3] for ps in mod_pumping_scheme]))
BF_keo=["" for i in range(sum_of_periods+1)]

for frame_ind in range(len(keo_fit_filenames)):
    fn=keo_fit_filenames[frame_ind].split("/")[-1]
    #
    
    ps_index, period_index, period_offset = get_position_on_scheme(keo_date_obs[frame_ind]+datetime.timedelta(seconds=keo_exp_sec[frame_ind]/2),mod_pumping_scheme)
    if ps_index==-1:
        BF_keo_ind=0
    else:
        BF_keo_ind=1
#         for i in range(ps_index):
        BF_keo_ind+=int(sum([ps[3] for ps in mod_pumping_scheme[0:ps_index]]))
        BF_keo_ind+=period_index
    #
    frame_date_axe, frame_x_dates = get_frame_time_axe(keo_date_obs[frame_ind], keo_exp_sec[frame_ind])
    overlap_coef=get_overlap_coef(frame_date_axe, mod_date_axe_clean, mod_clean_pumping)
    if overlap_coef<0.33 and period_offset<=1.:
        for i in range(len(keo_spcal_coefs)):
            if keo_spcal_coefs[i][0]==fn:
                i_find=i
                break
        if abs(keo_spcal_coefs[i_find][1]-keo_spcal_day_coef)<=3*keo_spcal_std:
            BF_keo[BF_keo_ind]=fn
            sign_flag=True
    i_find=-1
    sign_flag=False


# In[10]:

# BF_keo
save_filename='keo_140824_base.frames'
fid=open(save_filename,'w')
for bf in BF_keo:
    fid.write(bf+'\n')
fid.close()


# In[11]:

fig=plt.figure(figsize=(12.8,7.2))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
for i in range(len(clean_pumping)):
    plt.plot(x_dates_clean[i],clean_pumping[i],'r')
# for i in range(len(mod_clean_pumping)):
#     plt.plot(mod_x_dates_clean[i],mod_clean_pumping[i],'b')
plt.plot([],[],'r',lw=3,label='pump')
plt.plot([],[],'g',lw=3,label='sbig')
plt.plot([],[],'k',lw=3,label='keo')
for bf in BF_s1c:
    if bf!="":
        BF_s1c_date_obs=s1c_get_date_obs(s1c_fit_path+"/"+bf,ut_shift=-4)
        BF_s1c_exp_sec=s1c_get_exp_sec(s1c_fit_path+"/"+bf)
        BF_s1c_date_axe, BF_s1c_x_dates = get_frame_time_axe(BF_s1c_date_obs, BF_s1c_exp_sec)
        plt.plot([BF_s1c_x_dates[0], BF_s1c_x_dates[-1]],[0.7, 0.7],'g',lw=3)
for bf in BF_keo:
    if bf!="":
        BF_keo_date_obs=keo_get_date_obs(keo_fit_path+"/"+bf)
        BF_keo_exp_sec=keo_get_exp_sec(keo_fit_path+"/"+bf)
        BF_keo_date_axe, BF_keo_x_dates = get_frame_time_axe(BF_keo_date_obs, BF_keo_exp_sec)
        plt.plot([BF_keo_x_dates[0], BF_keo_x_dates[-1]],[0.8, 0.8],'k',lw=3)       
plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2014,8,24,17,50))
xlim_right=dates.date2num(datetime.datetime(2014,8,24,19,46))
# xlim_left=dates.date2num(datetime.datetime(2014,8,24,17,48))
# xlim_right=dates.date2num(datetime.datetime(2014,8,24,18,10))
ax.set_xlim(xlim_left,xlim_right)
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.title("14/08/24 Base frames' expositions")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.grid()
# plt.show()
plt.savefig('140824_base_frames.png')
plt.close()


# # 14/08/26

# In[12]:

# aug 26/08/14
pumping_scheme_file="pump140826.scheme"

s1c_fit_path="../data/140826/s1c"
keo_fit_path="../data/140826/keo"

s1c_spcal_fname = "../spectrophotometric_calibration/s1c_140826.spcal"
s1c_spcal_day_fname = "../spectrophotometric_calibration/s1c_140826_day.spcal"
s1c_spcal_day_coef, s1c_spcal_std=get_spcal_day_coefs(s1c_spcal_day_fname)    
s1c_spcal_coefs=get_spcal_coefs(s1c_spcal_fname)   

keo_spcal_fname = "../spectrophotometric_calibration/keo_140826.spcal"
keo_spcal_day_fname = "../spectrophotometric_calibration/keo_140826_day.spcal"
keo_spcal_day_coef, keo_spcal_std=get_spcal_day_coefs(keo_spcal_day_fname)
keo_spcal_coefs=get_spcal_coefs(keo_spcal_fname)

prohibited_area_min=2.
pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)
mod_pumping_scheme=get_mod_pumping_scheme(pumping_scheme,prohibited_area_min)


# In[13]:

s1c_fit_filenames=sorted([s1c_fit_path+'/'+fn for fn in next(os.walk(s1c_fit_path))[2]])
keo_fit_filenames=sorted([keo_fit_path+'/'+fn for fn in next(os.walk(keo_fit_path))[2]])

s1c_date_obs=[s1c_get_date_obs(f,ut_shift=-4) for f in s1c_fit_filenames]
s1c_exp_sec=[s1c_get_exp_sec(f) for f in s1c_fit_filenames]

keo_date_obs=[keo_get_date_obs(f) for f in keo_fit_filenames]
keo_exp_sec=[keo_get_exp_sec(f) for f in keo_fit_filenames]

date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
mod_date_axe_clean, mod_x_dates_clean, mod_clean_pumping = make_clean_pumping_scheme(mod_pumping_scheme)


# In[14]:

sum_of_periods=int(sum([ps[3] for ps in mod_pumping_scheme]))
BF_s1c=["" for i in range(sum_of_periods+1)]

for frame_ind in range(len(s1c_fit_filenames)):
    fn=s1c_fit_filenames[frame_ind].split("/")[-1]
    #
    
    ps_index, period_index, period_offset = get_position_on_scheme(s1c_date_obs[frame_ind]+datetime.timedelta(seconds=s1c_exp_sec[frame_ind]/2),mod_pumping_scheme)
    if ps_index==-1:
        BF_s1c_ind=0
    else:
        BF_s1c_ind=1
#         for i in range(ps_index-1):
        BF_s1c_ind+=int(sum([ps[3] for ps in mod_pumping_scheme[0:ps_index]]))
#         if frame_ind==689:
#                 print(BF_s1c_ind)
        BF_s1c_ind+=period_index
    #
    frame_date_axe, frame_x_dates = get_frame_time_axe(s1c_date_obs[frame_ind], s1c_exp_sec[frame_ind])
    overlap_coef=get_overlap_coef(frame_date_axe, mod_date_axe_clean, mod_clean_pumping)
    if overlap_coef<0.33 and period_offset<=1.:
        for i in range(len(s1c_spcal_coefs)):
            if s1c_spcal_coefs[i][0]==fn:
                i_find=i
                break
#                 if frame_ind==82: print("WoW!",abs(s1c_spcal_coefs[i][1]-s1c_spcal_day_coef),2*s1c_spcal_std);
        if abs(s1c_spcal_coefs[i_find][1]-s1c_spcal_day_coef)<=3*s1c_spcal_std:
#             print(frame_ind,BF_s1c_ind)
            BF_s1c[BF_s1c_ind]=fn
            sign_flag=True
#     if i_find>-1: print(frame_ind,fn,s1c_date_obs[frame_ind],overlap_coef,period_offset,sign_flag,abs(s1c_spcal_coefs[i_find][1]-s1c_spcal_day_coef),2*s1c_spcal_std);
    i_find=-1
    sign_flag=False


# In[15]:

# BF_s1c
save_filename='s1c_140826_base.frames'
fid=open(save_filename,'w')
for bf in BF_s1c:
    fid.write(bf+'\n')
fid.close()


# In[16]:

sum_of_periods=int(sum([ps[3] for ps in mod_pumping_scheme]))
BF_keo=["" for i in range(sum_of_periods+1)]

for frame_ind in range(len(keo_fit_filenames)):
    fn=keo_fit_filenames[frame_ind].split("/")[-1]
    #
    
    ps_index, period_index, period_offset = get_position_on_scheme(keo_date_obs[frame_ind]+datetime.timedelta(seconds=keo_exp_sec[frame_ind]/2),mod_pumping_scheme)
    if ps_index==-1:
        BF_keo_ind=0
    else:
        BF_keo_ind=1
#         for i in range(ps_index):
        BF_keo_ind+=int(sum([ps[3] for ps in mod_pumping_scheme[0:ps_index]]))
        BF_keo_ind+=period_index
    #
    frame_date_axe, frame_x_dates = get_frame_time_axe(keo_date_obs[frame_ind], keo_exp_sec[frame_ind])
    overlap_coef=get_overlap_coef(frame_date_axe, mod_date_axe_clean, mod_clean_pumping)
    if overlap_coef<0.33 and period_offset<=1.:
        for i in range(len(keo_spcal_coefs)):
            if keo_spcal_coefs[i][0]==fn:
                i_find=i
                break
        if abs(keo_spcal_coefs[i_find][1]-keo_spcal_day_coef)<=3*keo_spcal_std:
            BF_keo[BF_keo_ind]=fn
            sign_flag=True
    i_find=-1
    sign_flag=False


# In[17]:

# BF_keo
save_filename='keo_140826_base.frames'
fid=open(save_filename,'w')
for bf in BF_keo:
    fid.write(bf+'\n')
fid.close()


# In[18]:

fig=plt.figure(figsize=(12.8,7.2))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(2,62,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
for i in range(len(clean_pumping)):
    plt.plot(x_dates_clean[i],clean_pumping[i],'r')
# for i in range(len(mod_clean_pumping)):
#     plt.plot(mod_x_dates_clean[i],mod_clean_pumping[i],'b')
plt.plot([],[],'r',lw=3,label='pump')
plt.plot([],[],'g',lw=3,label='s1c')
plt.plot([],[],'k',lw=3,label='keo')
for bf in BF_s1c:
    if bf!="":        
        BF_s1c_date_obs=s1c_get_date_obs(s1c_fit_path+"/"+bf,ut_shift=-4)
        BF_s1c_exp_sec=s1c_get_exp_sec(s1c_fit_path+"/"+bf)        
        BF_s1c_date_axe, BF_s1c_x_dates = get_frame_time_axe(BF_s1c_date_obs, BF_s1c_exp_sec)
#         if bf==BF_s1c[-1]: print(BF_s1c_date_axe[0],BF_s1c_date_axe[-1]);
#         if bf==BF_s1c[-1]: print(BF_s1c_x_dates[0], BF_s1c_x_dates[-1]);

        plt.plot([BF_s1c_x_dates[0], BF_s1c_x_dates[-1]],[0.7, 0.7],'g',lw=3)
for bf in BF_keo:
    if bf!="":
        BF_keo_date_obs=keo_get_date_obs(keo_fit_path+"/"+bf)
        BF_keo_exp_sec=keo_get_exp_sec(keo_fit_path+"/"+bf)
        BF_keo_date_axe, BF_keo_x_dates = get_frame_time_axe(BF_keo_date_obs, BF_keo_exp_sec)
        plt.plot([BF_keo_x_dates[0], BF_keo_x_dates[-1]],[0.8, 0.8],'k',lw=3)       
plt.xticks(rotation=90)
plt.grid()
xlim_left=dates.date2num(datetime.datetime(2014,8,26,18,10))
xlim_right=dates.date2num(datetime.datetime(2014,8,26,21,20))
ax.set_xlim(xlim_left,xlim_right)
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.title("14/08/26 Base frames' expositions")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
# plt.show()
plt.savefig('140826_base_frames.png')
plt.close()


# # 16/08/29

# In[19]:

# aug 29/08/16
pumping_scheme_file="pump160829.scheme"

sbig_fit_path="../data/160829/sbig"
keo_fit_path="../data/160829/keo"

sbig_spcal_fname = "../spectrophotometric_calibration/sbig_160829.spcal"
sbig_spcal_day_fname = "../spectrophotometric_calibration/sbig_160829_day.spcal"
sbig_spcal_day_coef, sbig_spcal_std=get_spcal_day_coefs(sbig_spcal_day_fname)    
sbig_spcal_coefs=get_spcal_coefs(sbig_spcal_fname)   

keo_spcal_fname = "../spectrophotometric_calibration/keo_160829.spcal"
keo_spcal_day_fname = "../spectrophotometric_calibration/keo_160829_day.spcal"
keo_spcal_day_coef, keo_spcal_std=get_spcal_day_coefs(keo_spcal_day_fname)
keo_spcal_coefs=get_spcal_coefs(keo_spcal_fname)

prohibited_area_min=2.
pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)
mod_pumping_scheme=get_mod_pumping_scheme(pumping_scheme,prohibited_area_min)


# In[20]:

sbig_spcal_coefs_new=[]
for i in range(len(sbig_spcal_coefs)):
    if sbig_spcal_coefs[i][1]==0.0:
        sbig_spcal_coefs_new.append((sbig_spcal_coefs[i][0],sbig_spcal_day_coef))        
    else:
        sbig_spcal_coefs_new.append(sbig_spcal_coefs[i])
sbig_spcal_coefs=sbig_spcal_coefs_new
# sbig_spcal_coefs


# In[21]:

sbig_fit_filenames=sorted([sbig_fit_path+'/'+fn for fn in next(os.walk(sbig_fit_path))[2]])
keo_fit_filenames=sorted([keo_fit_path+'/'+fn for fn in next(os.walk(keo_fit_path))[2]])

sbig_date_obs=[sbig_get_date_obs(f) for f in sbig_fit_filenames]
sbig_exp_sec=[sbig_get_exp_sec(f) for f in sbig_fit_filenames]

keo_date_obs=[keo_get_date_obs(f) for f in keo_fit_filenames]
keo_exp_sec=[keo_get_exp_sec(f) for f in keo_fit_filenames]

date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
mod_date_axe_clean, mod_x_dates_clean, mod_clean_pumping = make_clean_pumping_scheme(mod_pumping_scheme)


# In[22]:

sum_of_periods=int(sum([ps[3] for ps in mod_pumping_scheme]))
BF_sbig=["" for i in range(sum_of_periods+1)]

for frame_ind in range(len(sbig_fit_filenames)):
    fn=sbig_fit_filenames[frame_ind].split("/")[-1]
    #
    
    ps_index, period_index, period_offset = get_position_on_scheme(sbig_date_obs[frame_ind]+datetime.timedelta(seconds=sbig_exp_sec[frame_ind]/2),mod_pumping_scheme)
    if ps_index==-1:
        BF_sbig_ind=0
    else:
        BF_sbig_ind=1
#         for i in range(ps_index-1):
        BF_sbig_ind+=int(sum([ps[3] for ps in mod_pumping_scheme[0:ps_index]]))
#         if frame_ind==689:
#                 print(BF_sbig_ind)
        BF_sbig_ind+=period_index
    #
    frame_date_axe, frame_x_dates = get_frame_time_axe(sbig_date_obs[frame_ind], sbig_exp_sec[frame_ind])
    overlap_coef=get_overlap_coef(frame_date_axe, mod_date_axe_clean, mod_clean_pumping)
    if overlap_coef<0.33 and period_offset<=1.:
        for i in range(len(sbig_spcal_coefs)):
            if sbig_spcal_coefs[i][0]==fn:
                i_find=i
                break
#                 if frame_ind==82: print("WoW!",abs(sbig_spcal_coefs[i][1]-sbig_spcal_day_coef),2*sbig_spcal_std);
        if abs(sbig_spcal_coefs[i_find][1]-sbig_spcal_day_coef)<=3*sbig_spcal_std:
#             print(frame_ind,BF_sbig_ind)
            BF_sbig[BF_sbig_ind]=fn
            sign_flag=True
#     if i_find>-1: print(frame_ind,fn,sbig_date_obs[frame_ind],overlap_coef,period_offset,sign_flag,abs(sbig_spcal_coefs[i_find][1]-sbig_spcal_day_coef),2*sbig_spcal_std);
    i_find=-1
    sign_flag=False


# In[23]:

# BF_sbig
save_filename='sbig_160829_base.frames'
fid=open(save_filename,'w')
for bf in BF_sbig:
    fid.write(bf+'\n')
fid.close()


# In[24]:

sum_of_periods=int(sum([ps[3] for ps in mod_pumping_scheme]))
BF_keo=["" for i in range(sum_of_periods+1)]

for frame_ind in range(len(keo_fit_filenames)):
    fn=keo_fit_filenames[frame_ind].split("/")[-1]
    #
    
    ps_index, period_index, period_offset = get_position_on_scheme(keo_date_obs[frame_ind]+datetime.timedelta(seconds=keo_exp_sec[frame_ind]/2),mod_pumping_scheme)
    if ps_index==-1:
        BF_keo_ind=0
    else:
        BF_keo_ind=1
#         for i in range(ps_index):
        BF_keo_ind+=int(sum([ps[3] for ps in mod_pumping_scheme[0:ps_index]]))
        BF_keo_ind+=period_index
    #
    frame_date_axe, frame_x_dates = get_frame_time_axe(keo_date_obs[frame_ind], keo_exp_sec[frame_ind])
    overlap_coef=get_overlap_coef(frame_date_axe, mod_date_axe_clean, mod_clean_pumping)
    if overlap_coef<0.33 and period_offset<=1.:
        for i in range(len(keo_spcal_coefs)):
            if keo_spcal_coefs[i][0]==fn:
                i_find=i
                break
        if abs(keo_spcal_coefs[i_find][1]-keo_spcal_day_coef)<=3*keo_spcal_std:
            BF_keo[BF_keo_ind]=fn
            sign_flag=True
    i_find=-1
    sign_flag=False


# In[25]:

# BF_keo
save_filename='keo_160829_base.frames'
fid=open(save_filename,'w')
for bf in BF_keo:
    fid.write(bf+'\n')
fid.close()


# In[26]:

fig=plt.figure(figsize=(12.8,7.2))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
# ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
for i in range(len(clean_pumping)):
    plt.plot(x_dates_clean[i],clean_pumping[i],'r')
# for i in range(len(mod_clean_pumping)):
#     plt.plot(mod_x_dates_clean[i],mod_clean_pumping[i],'b')
plt.plot([],[],'r',lw=3,label='pump')
plt.plot([],[],'g',lw=3,label='sbig')
plt.plot([],[],'k',lw=3,label='keo')
for bf in BF_sbig:
    if bf!="":        
        BF_sbig_date_obs=sbig_get_date_obs(sbig_fit_path+"/"+bf)
        BF_sbig_exp_sec=sbig_get_exp_sec(sbig_fit_path+"/"+bf)        
        BF_sbig_date_axe, BF_sbig_x_dates = get_frame_time_axe(BF_sbig_date_obs, BF_sbig_exp_sec)
        plt.plot([BF_sbig_x_dates[0], BF_sbig_x_dates[-1]],[0.7, 0.7],'g',lw=3)
for bf in BF_keo:
    if bf!="":
        BF_keo_date_obs=keo_get_date_obs(keo_fit_path+"/"+bf)
        BF_keo_exp_sec=keo_get_exp_sec(keo_fit_path+"/"+bf)
        BF_keo_date_axe, BF_keo_x_dates = get_frame_time_axe(BF_keo_date_obs, BF_keo_exp_sec)
        plt.plot([BF_keo_x_dates[0], BF_keo_x_dates[-1]],[0.8, 0.8],'k',lw=3)   
plt.grid()
plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2016,8,29,17,28,0))
xlim_right=dates.date2num(datetime.datetime(2016,8,29,20,25,0))
ax.set_xlim(xlim_left,xlim_right)
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.title("16/08/29 Base frames' expositions")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
# plt.show()
plt.savefig('160829_base_frames.png')
plt.close()

