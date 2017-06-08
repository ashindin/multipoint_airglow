
# coding: utf-8

# In[1]:

import datetime, os
import numpy as np
import scipy.io as si


# In[2]:

from matplotlib import dates
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# In[3]:

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

def get_photom_start_time(par_filename,ut_shift=-3):
    fid=open(par_filename,'r', encoding='cp1251')
    s=fid.readline()[-21:-1]
    fid.close()    
    bin_filename=par_filename[0:-4]+".dat"
    num_sec = int(os.stat(bin_filename).st_size/2/3/1000)    
    d=datetime.datetime.strptime(s, "%b %d %H:%M:%S %Y")+datetime.timedelta(hours=ut_shift, seconds=-num_sec)
    return d

def get_photom_data(filename):
    data=np.fromfile(filename, dtype=np.int16).reshape((-1,3))
    return data[:,2]

def get_photom_time_axe(data_len,photom_start_date):
    num_sec=int(data_len/1000)
    date_axe=[photom_start_date+datetime.timedelta(seconds=i) for i in range(num_sec)]
    x_dates=dates.date2num(date_axe)
    return date_axe, x_dates

def get_hf_data(filename,st_date):
    data_dic=si.loadmat(filenames[0])
    data=data_dic['HF_Data_1sek'][0]
    sec_axe=data_dic['HF_UTCSecScale_1sek']
    num_sec=len(data)
    date_axe=[st_date+datetime.timedelta(seconds=int(sec_axe[0][i])) for i in range(num_sec)]
    x_dates=dates.date2num(date_axe)    
    return data, date_axe, x_dates

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


# In[4]:

# aug 24/08/14
ut_shift=-4
filenames=[
"../data/140824/photometers/24081401.dat",
"../data/140824/photometers/24081402.dat",
"../data/140824/photometers/24081403.dat",
"../data/140824/photometers/24081404.dat",
]
pumping_scheme_file="pump140824.scheme"


# In[5]:

pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)

par_filenames=[f[0:-4]+".par" for f in filenames]
photom_start_dates=[get_photom_start_time(par_filename,-4) for par_filename in par_filenames]

date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
fig=plt.figure(figsize=(16,6))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
for ind in range(len(filenames)):
    data=get_photom_data(filenames[ind])
    date_axe, x_dates = get_photom_time_axe(len(data),photom_start_dates[ind])
    plt.plot(x_dates,2*data[::1000]/np.max(data),'b')
for i in range(len(clean_pumping)):
    plt.plot(x_dates_clean[i],clean_pumping[i],'r')
plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2014,8,24,17,50))
xlim_right=dates.date2num(datetime.datetime(2014,8,24,19,45))
ax.set_xlim(xlim_left,xlim_right)
plt.grid()
plt.title("exam_pumping_scheme_140824.png")
plt.savefig("exam_pumping_scheme_140824.png")
plt.close()
# plt.show()


# In[6]:

# aug 26/08/14
ut_shift=-4
filenames=[
"../data/140826/photometers/26081401.dat",
"../data/140826/photometers/26081402.dat",
"../data/140826/photometers/26081403.dat",
"../data/140826/photometers/26081404.dat",
"../data/140826/photometers/26081405.dat",
"../data/140826/photometers/26081406.dat",
]
pumping_scheme_file="pump140826.scheme"


# In[7]:

pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)

par_filenames=[f[0:-4]+".par" for f in filenames]
photom_start_dates=[get_photom_start_time(par_filename,-4) for par_filename in par_filenames]

date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
fig=plt.figure(figsize=(16,6))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(2,62,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
for ind in range(len(filenames)):
    data=get_photom_data(filenames[ind])
    date_axe, x_dates = get_photom_time_axe(len(data),photom_start_dates[ind])
    plt.plot(x_dates,2*data[::1000]/np.max(data),'b')
for i in range(len(clean_pumping)):
    plt.plot(x_dates_clean[i],clean_pumping[i],'r')
plt.grid()
plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2014,8,26,18,10))
xlim_right=dates.date2num(datetime.datetime(2014,8,26,21,15))
ax.set_xlim(xlim_left,xlim_right)
plt.title("exam_pumping_scheme_140826.png")
plt.savefig("exam_pumping_scheme_140826.png")
plt.close()
# plt.show()


# In[8]:

# aug 29/08/16
filenames=["../data/160829/hf/HF_Data_OneSec0829.mat"]
st_date=datetime.datetime(2016,8,29)
pumping_scheme_file="pump160829.scheme"


# In[9]:

pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)

date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
fig=plt.figure(figsize=(16,6))
ax=plt.axes()
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
# ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)
for ind in range(len(filenames)):
    data, date_axe, x_dates = get_hf_data(filenames[ind],st_date)
    plt.plot(x_dates, data/np.max(data),'b')
for i in range(len(clean_pumping)):
    plt.plot(x_dates_clean[i],0.35*clean_pumping[i],'r')
plt.grid()
plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2016,8,29,17,0,0))
xlim_right=dates.date2num(datetime.datetime(2016,8,29,20,30,0))
ax.set_xlim(xlim_left,xlim_right)
ax.set_ylim(0.1,0.5)
plt.title("exam_pumping_scheme_160829.png")
plt.savefig("exam_pumping_scheme_160829.png")
plt.close()
# plt.show()

