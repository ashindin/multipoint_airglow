
# coding: utf-8

# In[1]:

import datetime, os, sys, inspect
import numpy as np


# In[2]:

from matplotlib import dates
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


coef=5.15/5.627;

# In[3]:

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../comparing")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)


# In[4]:

from get_ranges import *


# In[5]:

replace_str='scp ip2165s:~/git/multipoint_airglow'


# In[6]:

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

def get_model_height(dat_filenames):
    hei_xdates=[]
    hei=[]
    for i in range(len(dat_filenames)):
        dat_fn=dat_filenames[i]
        fid=open(dat_fn,'r')
        fid.readline()
        success_str=fid.readline()
        status_str=fid.readline()
        if success_str!="SUCCESS = True\n" or status_str!="STATUS = 0\n":
            fid.close()
            continue
        fid.readline()
        fid.readline()
        fid.readline()
        hei.append(float(fid.readline().split(' ')[2]))

        date_str=dat_fn.split('/')[-1][0:20]
        yr=int(date_str[0:2])+2000
        mon=int(date_str[2:4])
        day=int(date_str[4:6])
        hh=int(date_str[7:9])
        mm=int(date_str[9:11])
        ss=int(date_str[11:13])
        mks=int(date_str[14:20])
        # print(date_str)
        # print(yr,mon,day,hh,mm,ss,mks)
        hei_xdates.append(dates.date2num(datetime.datetime(yr,mon,day,hh,mm,ss,mks)))
    #     print(xdate, hei)
        fid.close()
    return hei_xdates, hei

def get_model_height2(dat_filenames,smatrix_bool,stype):
    
    dat_filenames2=[]
    for i in range(len(dat_filenames)):
        dat_fn=dat_filenames[i].replace('sphere',stype)
        dat_filenames2.append(dat_fn)
    dat_filenames=dat_filenames2
    
#     print(dat_filenames)
    
    hei_xdates=np.zeros(len(dat_filenames))
    hei=np.zeros(len(dat_filenames))
    
    if stype=='sphere': sval=0;
    if stype=='spheroid': sval=1;
    if stype=='ellipsoid': sval=2;
    if stype=='spheroid_incl': sval=3;
    if stype=='drop': sval=4;
    
    for i in range(len(dat_filenames)):
        
        if smatrix_bool[i,sval]==False:
            hei[i]=np.nan
            hei_xdates[i]=np.nan
            continue
        
        dat_fn=dat_filenames[i]
        if os.path.isfile(dat_fn)==False:
            hei[i]=np.nan
            hei_xdates[i]=np.nan
            continue
        
        fid=open(dat_fn,'r')
        fid.readline()
        success_str=fid.readline()
        status_str=fid.readline()
        if success_str!="SUCCESS = True\n" or status_str!="STATUS = 0\n":
            fid.close()
            hei[i]=np.nan
            hei_xdates[i]=np.nan
            continue
        fid.readline()
        fid.readline()
        fid.readline()
        hei[i]=float(fid.readline().split(' ')[2])

        date_str=dat_fn.split('/')[-1][0:20]
        yr=int(date_str[0:2])+2000
        mon=int(date_str[2:4])
        day=int(date_str[4:6])
        hh=int(date_str[7:9])
        mm=int(date_str[9:11])
        ss=int(date_str[11:13])
        mks=int(date_str[14:20])
        # print(date_str)
        # print(yr,mon,day,hh,mm,ss,mks)
        hei_xdates[i]=dates.date2num(datetime.datetime(yr,mon,day,hh,mm,ss,mks))
    #     print(xdate, hei)
        fid.close()
    return hei_xdates, hei

def get_model_height3(dat_filenames,smatrix_bool,stype, x_dates_clean, clean_pumping):
    
    dat_filenames2=[]
    for i in range(len(dat_filenames)):
        dat_fn=dat_filenames[i].replace('sphere',stype)
        dat_filenames2.append(dat_fn)
    dat_filenames=dat_filenames2
    
#     print(dat_filenames)
    
    hei_xdates=np.zeros(len(dat_filenames))
    hei=np.zeros(len(dat_filenames))
    dzb=np.zeros(len(dat_filenames))
    dzt=np.zeros(len(dat_filenames))
    
    if stype=='sphere': sval=0;
    if stype=='spheroid': sval=1;
    if stype=='ellipsoid': sval=2;
    if stype=='spheroid_incl': sval=3;
    if stype=='drop': sval=4;
    
    for i in range(len(dat_filenames)):
        
        if smatrix_bool[i,sval]==False:
            hei[i]=np.nan
            hei_xdates[i]=np.nan
            continue
        
        dat_fn=dat_filenames[i]
        if os.path.isfile(dat_fn)==False:
            hei[i]=np.nan
            hei_xdates[i]=np.nan
            continue
        
        fid=open(dat_fn,'r')
        fid.readline()
        success_str=fid.readline()
        status_str=fid.readline()
        if success_str!="SUCCESS = True\n" or status_str!="STATUS = 0\n":
            fid.close()
            hei[i]=np.nan
            dzb[i]=np.nan
            dzt[i]=np.nan
            hei_xdates[i]=np.nan
            continue
        fid.readline()
        fid.readline()
        fid.readline()
        hei[i]=float(fid.readline().split(' ')[2])

        date_str=dat_fn.split('/')[-1][0:20]
        yr=int(date_str[0:2])+2000
        mon=int(date_str[2:4])
        day=int(date_str[4:6])
        hh=int(date_str[7:9])
        mm=int(date_str[9:11])
        ss=int(date_str[11:13])
        mks=int(date_str[14:20])
        # print(date_str)
        # print(yr,mon,day,hh,mm,ss,mks)
        hei_xdates[i]=dates.date2num(datetime.datetime(yr,mon,day,hh,mm,ss,mks))
        
#         print(hei_xdates[i])
#         print(x_dates_clean, clean_pumping)
        val=0.
        for j in range(len(clean_pumping)):
#             if 
            temp_val=np.interp(hei_xdates[i], x_dates_clean[j], clean_pumping[j])
#             print(temp_val)
            if np.isnan(temp_val)==False:
                val+=temp_val
#                 print('wow')
#         print(val)
        if val==0.:
            hei[i]=np.nan
            hei_xdates[i]=np.nan
        
        fid.close()
        if sval==3:
            dzb[i], dzt[i] = get_spheroid_incl_range(dat_fn)
        elif sval==2:
            dzb[i], dzt[i] = get_ellipsoid_range(dat_fn)
        elif sval==4:
#             print(dat_fn)
            dzb[i], dzt[i] = get_drop_range(dat_fn)
        else:
            dzb[i]=0.
            dzt[i]=0.
        
    return hei_xdates, hei, dzb, dzt

def get_model_fun(dat_fn):
    if os.path.isfile(dat_fn)==False:
        return np.nan
    fid=open(dat_fn,'r')
    fid.readline()
    success_str=fid.readline()
    status_str=fid.readline()
    if success_str!="SUCCESS = True\n" or status_str!="STATUS = 0\n":
        fid.close()
        return np.nan
    fid.readline()
    fid.readline()
    fun=float(fid.readline().split('=')[1])
    return fun
# get_model_fun('d:/git/multipoint_airglow/solving/drop/140824_175227_000000_drop.dat')


# In[7]:

def get_model_pars(dat_fn):
    if os.path.isfile(dat_fn)==False:
        return np.nan
    fid=open(dat_fn,'r')
    fid.readline()
    success_str=fid.readline()
    status_str=fid.readline()
#     print(status_str)
    if success_str!="SUCCESS = True\n" or status_str!="STATUS = 0\n":
        fid.close()
        return np.nan
    fid.readline()
    fid.readline()
    fun=float(fid.readline().split('=')[1])
    pars_list=fid.readline().split()
    pars=[float(pl) for pl in pars_list]
    mod_pos=pars[:3]
    mod_pos[0]*=np.pi/180
    mod_pos[1]*=np.pi/180
    return fun, mod_pos, pars[3::]


# In[8]:

day_str='140824'
pumping_scheme_file="../image_processing/pump140824.scheme"
ion_fn='../comparing/refl_heights_140824.table'


# In[9]:

pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)
date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
pumping_scheme


# In[10]:

dat_sphere_filenames=[]
for fn in next(os.walk('../solving/sphere/'))[2]:
    if fn[-3:]=='dat' and fn[0:6]==day_str:
        dat_sphere_filenames.append('../solving/sphere/'+fn)
dat_sphere_filenames=sorted(dat_sphere_filenames)


# In[11]:

smatrix=np.zeros((len(dat_sphere_filenames),5))
smatrix_bool=np.zeros((len(dat_sphere_filenames),5),dtype=np.bool)
for i in range(len(dat_sphere_filenames)):
    dat_fn=dat_sphere_filenames[i]
#     print(dat_fn,dat_fn.replace('sphere','spheroid'))
    smatrix[i,0]=get_model_fun(dat_fn)
    smatrix[i,1]=get_model_fun(dat_fn.replace('sphere','spheroid'))
    smatrix[i,2]=get_model_fun(dat_fn.replace('sphere','ellipsoid'))
    smatrix[i,3]=get_model_fun(dat_fn.replace('sphere','spheroid_incl'))
    smatrix[i,4]=get_model_fun(dat_fn.replace('sphere','drop'))
    
    argmin=np.argmin(smatrix[i,:])
    smatrix_bool[i,argmin]=True
# smatrix_bool


# In[12]:

hei_sphere_xdates, hei_sphere, dzb_sphere, dzt_sphere = get_model_height3(dat_sphere_filenames,smatrix_bool,'sphere', x_dates_clean, clean_pumping)
hei_spheroid_xdates, hei_spheroid, dzb_spheroid, dzt_spheroid = get_model_height3(dat_sphere_filenames,smatrix_bool,'spheroid', x_dates_clean, clean_pumping)
hei_ellipsoid_xdates, hei_ellipsoid, dzb_ellipsoid, dzt_ellipsoid = get_model_height3(dat_sphere_filenames,smatrix_bool,'ellipsoid', x_dates_clean, clean_pumping)
hei_spheroid_incl_xdates, hei_spheroid_incl, dzb_spheroid_incl, dzt_spheroid_incl = get_model_height3(dat_sphere_filenames,smatrix_bool,'spheroid_incl', x_dates_clean, clean_pumping)
hei_drop_xdates, hei_drop, dzb_drop, dzt_drop = get_model_height3(dat_sphere_filenames,smatrix_bool,'drop', x_dates_clean, clean_pumping)


# In[13]:

fid_ion=open(ion_fn,'r')
ion_lines=fid_ion.readlines()[1::]
fid_ion.close()

ion_xdates=[]
ion_refl=[]
ion_b=[]
ion_t=[]
for i in range(len(ion_lines)):
    ion_line=ion_lines[i].split('\n')[0].split(',')
    ion_xdates.append(dates.date2num(datetime.datetime.strptime(ion_line[0],'%Y-%m-%dT%H:%M:%S')))
    ion_refl.append(float(ion_line[1])*1000)
    ion_b.append(float(ion_line[2])*1000)
    ion_t.append(float(ion_line[3])*1000)


# In[14]:

err_range=231

fig=plt.figure(figsize=(6.69,3.34))
ax=plt.axes(position=[0.1,0.22,0.87,0.75])
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)

for i in range(len(clean_pumping)):
    plt.plot(x_dates_clean[i],0+100*clean_pumping[i],'k')


conc=[]
for i in range(len(hei_ellipsoid)):
    if dzt_ellipsoid[i]+dzb_ellipsoid[i]>5000 and dzt_ellipsoid[i]+dzb_ellipsoid[i] < 150000 and hei_ellipsoid[i]>175000:
        dat_fn=dat_sphere_filenames[i].replace('sphere','ellipsoid')

        fun, mod_pos, pars = get_model_pars(dat_fn)
        conc.append(pars[0])
        plt.plot([hei_ellipsoid_xdates[i],hei_ellipsoid_xdates[i]],[coef*pars[0]-err_range,coef*pars[0]+err_range],'m',lw=1)
        plt.plot(hei_ellipsoid_xdates[i], coef*pars[0],'ms',mec='k')

for i in range(len(hei_spheroid_incl)):
    if dzt_spheroid_incl[i]+dzb_spheroid_incl[i]>5000 and dzt_spheroid_incl[i]+dzb_spheroid_incl[i] < 150000 and hei_spheroid_incl[i]>175000:
        dat_fn=dat_sphere_filenames[i].replace('sphere','spheroid_incl')
        fun, mod_pos, pars = get_model_pars(dat_fn)
        conc.append(pars[0])
        plt.plot([hei_spheroid_incl_xdates[i],hei_spheroid_incl_xdates[i]],[coef*pars[0]-err_range,coef*pars[0]+err_range],'r',lw=1)
        plt.plot(hei_spheroid_incl_xdates[i], coef*pars[0],'ro',mec='k')

for i in range(len(hei_drop)):
    if dzt_drop[i]+dzb_drop[i]>5000 and dzt_drop[i]+dzb_drop[i]<150000 and hei_drop[i]>175000:

        dat_fn=dat_sphere_filenames[i].replace('sphere','drop')
        fun, mod_pos, pars = get_model_pars(dat_fn)
        conc.append(pars[0])
        plt.plot([hei_drop_xdates[i],hei_drop_xdates[i]],[coef*pars[0]-err_range,coef*pars[0]+err_range],'g',lw=1)
        plt.plot(hei_drop_xdates[i], coef*pars[0],'gv',mec='k')

plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2014,8,24,17,50))
xlim_right=dates.date2num(datetime.datetime(2014,8,24,19,45))
ax.set_xlim(xlim_left,xlim_right)
ax.set_ylim(0,2500)
plt.xlabel("Время ЧЧ:ММ, UTC")
plt.ylabel("$p_1$, см$^{-3}$")
plt.grid()


h_sp, = plt.plot(0,0,'ro',mec='k',label='M4') # spheroid_incl
h_dr, = plt.plot(0,0,'gv',mec='k',label='M5') # drop
h_el, = plt.plot(0,0,'ms',mec='k',label='M3') # ellipsoid
plt.legend([h_el,h_sp,h_dr],['M3','M4','M5'])


# plt.title("Airglow's center height 140824")
plt.savefig("fig9_1.png",dpi=300)
plt.savefig("fig9_1.pdf",dpi=300)
plt.savefig("fig9_1.eps",dpi=300)
plt.close()
# plt.show()

print(np.mean(coef*np.array(conc)),np.std(coef*np.array(conc)))


# In[15]:

np.cos(12*np.pi/180)*4300


# In[16]:

day_str='140826'
pumping_scheme_file="../image_processing/pump140826.scheme"
ion_fn='../comparing/refl_heights_140826.table'


# In[17]:

pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)
date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
pumping_scheme


# In[18]:

dat_sphere_filenames=[]
for fn in next(os.walk('../solving/sphere/'))[2]:
    if fn[-3:]=='dat' and fn[0:6]==day_str:
        dat_sphere_filenames.append('../solving/sphere/'+fn)
dat_sphere_filenames=sorted(dat_sphere_filenames)


# In[19]:

smatrix=np.zeros((len(dat_sphere_filenames),5))
smatrix_bool=np.zeros((len(dat_sphere_filenames),5),dtype=np.bool)
for i in range(len(dat_sphere_filenames)):
    dat_fn=dat_sphere_filenames[i]
#     print(dat_fn,dat_fn.replace('sphere','spheroid'))
    smatrix[i,0]=get_model_fun(dat_fn)
    smatrix[i,1]=get_model_fun(dat_fn.replace('sphere','spheroid'))
    smatrix[i,2]=get_model_fun(dat_fn.replace('sphere','ellipsoid'))
    smatrix[i,3]=get_model_fun(dat_fn.replace('sphere','spheroid_incl'))
    smatrix[i,4]=get_model_fun(dat_fn.replace('sphere','drop'))
    
    argmin=np.argmin(smatrix[i,:])
    smatrix_bool[i,argmin]=True
# smatrix_bool


# In[20]:

hei_sphere_xdates, hei_sphere, dzb_sphere, dzt_sphere  = get_model_height3(dat_sphere_filenames,smatrix_bool,'sphere', x_dates_clean, clean_pumping)
hei_spheroid_xdates, hei_spheroid, dzb_spheroid, dzt_spheroid  = get_model_height3(dat_sphere_filenames,smatrix_bool,'spheroid', x_dates_clean, clean_pumping)
hei_ellipsoid_xdates, hei_ellipsoid, dzb_ellipsoid, dzt_ellipsoid  = get_model_height3(dat_sphere_filenames,smatrix_bool,'ellipsoid', x_dates_clean, clean_pumping)
hei_spheroid_incl_xdates, hei_spheroid_incl, dzb_spheroid_incl, dzt_spheroid_incl  = get_model_height3(dat_sphere_filenames,smatrix_bool,'spheroid_incl', x_dates_clean, clean_pumping)
hei_drop_xdates, hei_drop, dzb_drop, dzt_drop  = get_model_height3(dat_sphere_filenames,smatrix_bool,'drop', x_dates_clean, clean_pumping)


# In[21]:

fid_ion=open(ion_fn,'r')
ion_lines=fid_ion.readlines()[1::]
fid_ion.close()

ion_xdates=[]
ion_refl=[]
ion_b=[]
ion_t=[]
for i in range(len(ion_lines)):
    ion_line=ion_lines[i].split('\n')[0].split(',')
    ion_xdates.append(dates.date2num(datetime.datetime.strptime(ion_line[0],'%Y-%m-%dT%H:%M:%S')))
    ion_refl.append(float(ion_line[1])*1000)
    ion_b.append(float(ion_line[2])*1000)
    ion_t.append(float(ion_line[3])*1000)


# In[22]:

err_range=231

fig=plt.figure(figsize=(6.69,3.34))
ax=plt.axes(position=[0.1,0.22,0.87,0.75])
ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(2,62,3)))
ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
xfmt = dates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(xfmt)

for i in range(len(clean_pumping)):
    plt.plot(x_dates_clean[i],0+100*clean_pumping[i],'k')

conc=[]
for i in range(len(hei_ellipsoid)):
    if dzt_ellipsoid[i]+dzb_ellipsoid[i]>5000 and dzt_ellipsoid[i]+dzb_ellipsoid[i] < 150000 and hei_ellipsoid[i]>175000:
        dat_fn=dat_sphere_filenames[i].replace('sphere','ellipsoid')
        fun, mod_pos, pars = get_model_pars(dat_fn)
        conc.append(pars[0])
        plt.plot([hei_ellipsoid_xdates[i],hei_ellipsoid_xdates[i]],[coef*pars[0]-err_range,coef*pars[0]+err_range],'g')
        plt.plot(hei_ellipsoid_xdates[i], coef*pars[0],'gv')

for i in range(len(hei_spheroid_incl)):
    if dzt_spheroid_incl[i]+dzb_spheroid_incl[i]>5000 and dzt_spheroid_incl[i]+dzb_spheroid_incl[i] < 150000 and hei_spheroid_incl[i]>175000 and hei_spheroid_incl[i]<310000:
        dat_fn=dat_sphere_filenames[i].replace('sphere','spheroid_incl')
        fun, mod_pos, pars = get_model_pars(dat_fn)
        conc.append(pars[0])
        plt.plot([hei_spheroid_incl_xdates[i],hei_spheroid_incl_xdates[i]],[coef*pars[0]-err_range,coef*pars[0]+err_range],'r',lw=1)
        plt.plot(hei_spheroid_incl_xdates[i], coef*pars[0],'ro',mec='k')

for i in range(len(hei_drop)):
    if dzt_drop[i]+dzb_drop[i]>5000 and dzt_drop[i]+dzb_drop[i]<150000 and hei_drop[i]>175000:
        dat_fn=dat_sphere_filenames[i].replace('sphere','drop')
        fun, mod_pos, pars = get_model_pars(dat_fn)
        conc.append(pars[0])
        plt.plot([hei_drop_xdates[i],hei_drop_xdates[i]],[coef*pars[0]-err_range,coef*pars[0]+err_range],'k')
        plt.plot(hei_drop_xdates[i], coef*pars[0],'k8')

plt.xticks(rotation=90)
xlim_left=dates.date2num(datetime.datetime(2014,8,26,19,26))
xlim_right=dates.date2num(datetime.datetime(2014,8,26,21,15))
ax.set_xlim(xlim_left,xlim_right)
ax.set_ylim(0,2500)
plt.xlabel("Время ЧЧ:ММ, UTC")
plt.ylabel("$p_1$, см$^{-3}$")
plt.grid()

h_sp, = plt.plot(0,0,'ro',mec='k',label='M4') # spheroid_incl
plt.legend([h_sp],['M4'])

# plt.title("Airglow's center height 140826")
# plt.savefig("fig9_2.png",dpi=300)
# plt.close()
# plt.show()

print(np.mean(coef*np.array(conc)),np.std(coef*np.array(conc)))


# In[23]:

# day_str='160829'
# pumping_scheme_file="../image_processing/pump160829.scheme"
# ion_fn='refl_heights_160829.table'


# In[24]:

# pumping_scheme=get_pumping_scheme_from_file(pumping_scheme_file)
# date_axe_clean, x_dates_clean, clean_pumping = make_clean_pumping_scheme(pumping_scheme)
# pumping_scheme


# In[25]:

# dat_sphere_filenames=[]
# for fn in next(os.walk('../solving/sphere/'))[2]:
#     if fn[-3:]=='dat' and fn[0:6]==day_str:
#         dat_sphere_filenames.append('../solving/sphere/'+fn)
# dat_sphere_filenames=sorted(dat_sphere_filenames)


# In[26]:

# smatrix=np.zeros((len(dat_sphere_filenames),5))
# smatrix_bool=np.zeros((len(dat_sphere_filenames),5),dtype=np.bool)
# for i in range(len(dat_sphere_filenames)):
#     dat_fn=dat_sphere_filenames[i]
# #     print(dat_fn,dat_fn.replace('sphere','spheroid'))
#     smatrix[i,0]=get_model_fun(dat_fn)
#     smatrix[i,1]=get_model_fun(dat_fn.replace('sphere','spheroid'))
#     smatrix[i,2]=get_model_fun(dat_fn.replace('sphere','ellipsoid'))
#     smatrix[i,3]=get_model_fun(dat_fn.replace('sphere','spheroid_incl'))
#     smatrix[i,4]=get_model_fun(dat_fn.replace('sphere','drop'))
    
#     argmin=np.argmin(smatrix[i,:])
#     smatrix_bool[i,argmin]=True
# # smatrix_bool


# In[27]:

# hei_sphere_xdates, hei_sphere, dzb_sphere, dzt_sphere = get_model_height3(dat_sphere_filenames,smatrix_bool,'sphere', x_dates_clean, clean_pumping)
# hei_spheroid_xdates, hei_spheroid, dzb_spheroid, dzt_spheroid = get_model_height3(dat_sphere_filenames,smatrix_bool,'spheroid', x_dates_clean, clean_pumping)
# hei_ellipsoid_xdates, hei_ellipsoid, dzb_ellipsoid, dzt_ellipsoid = get_model_height3(dat_sphere_filenames,smatrix_bool,'ellipsoid', x_dates_clean, clean_pumping)
# hei_spheroid_incl_xdates, hei_spheroid_incl, dzb_spheroid_incl, dzt_spheroid_incl = get_model_height3(dat_sphere_filenames,smatrix_bool,'spheroid_incl', x_dates_clean, clean_pumping)
# hei_drop_xdates, hei_drop, dzb_drop, dzt_drop = get_model_height3(dat_sphere_filenames,smatrix_bool,'drop', x_dates_clean, clean_pumping)


# In[28]:

# fid_ion=open(ion_fn,'r')
# ion_lines=fid_ion.readlines()[1::]
# fid_ion.close()

# ion_xdates=[]
# ion_refl=[]
# ion_uh=[]
# # ion_t=[]
# for i in range(len(ion_lines)):
#     ion_line=ion_lines[i].split('\n')[0].split(',')
#     ion_xdates.append(dates.date2num(datetime.datetime.strptime(ion_line[0],'%Y-%m-%dT%H:%M:%S')))
#     ion_refl.append(float(ion_line[1])*1000)
#     ion_uh.append(float(ion_line[2])*1000)
# #     ion_t.append(float(ion_line[3])*1000)


# In[29]:

# fig=plt.figure(figsize=(16,6))
# ax=plt.axes()
# ax.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
# # ax.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
# xfmt = dates.DateFormatter('%H:%M')
# ax.xaxis.set_major_formatter(xfmt)

# for i in range(len(clean_pumping)):
#     plt.plot(x_dates_clean[i],150000+50000*clean_pumping[i],'k')

# plt.plot(ion_xdates,ion_refl,'r-')
# plt.plot(ion_xdates,ion_uh,'g-')
# # plt.plot(ion_xdates,ion_t,'g-')

# plt.plot(hei_sphere_xdates, hei_sphere,'ro')
# plt.plot(hei_spheroid_xdates, hei_spheroid,'bs')

# for i in range(len(hei_ellipsoid)):
#     if dzt_ellipsoid[i]+dzb_ellipsoid[i]>5000 and dzt_ellipsoid[i]+dzb_ellipsoid[i] < 150000 and hei_ellipsoid[i]>175000 and i>131:
#         dat_fn=dat_sphere_filenames[i].replace('sphere','ellipsoid')
#         print(dat_fn.replace('dat','png').replace('..',replace_str) + ' .')
#         plt.plot([hei_ellipsoid_xdates[i],hei_ellipsoid_xdates[i]],[hei_ellipsoid[i]-dzb_ellipsoid[i],hei_ellipsoid[i]+dzt_ellipsoid[i]],'g')
#         plt.plot(hei_ellipsoid_xdates[i], hei_ellipsoid[i],'gv')

# for i in range(len(hei_spheroid_incl)):
#     if dzt_spheroid_incl[i]+dzb_spheroid_incl[i]>5000 and dzt_spheroid_incl[i]+dzb_spheroid_incl[i] < 150000 and hei_spheroid_incl[i]>175000 and i>131:
#         dat_fn=dat_sphere_filenames[i].replace('sphere','spheroid_incl')
#         print(dat_fn.replace('dat','png').replace('..',replace_str) + ' .')
#         plt.plot([hei_spheroid_incl_xdates[i],hei_spheroid_incl_xdates[i]],[hei_spheroid_incl[i]-dzb_spheroid_incl[i],hei_spheroid_incl[i]+dzt_spheroid_incl[i]],'m')
#         plt.plot(hei_spheroid_incl_xdates[i], hei_spheroid_incl[i],'mp')

# for i in range(len(hei_drop)):
#     if dzt_drop[i]+dzb_drop[i]>5000 and dzt_drop[i]+dzb_drop[i]<150000 and hei_drop[i]>175000 and i>131:
#         dat_fn=dat_sphere_filenames[i].replace('sphere','drop')
#         print(dat_fn.replace('dat','png').replace('..',replace_str) + ' .')
#         plt.plot([hei_drop_xdates[i],hei_drop_xdates[i]],[hei_drop[i]-dzb_drop[i],hei_drop[i]+dzt_drop[i]],'k')
#         plt.plot(hei_drop_xdates[i], hei_drop[i],'k8')

# plt.xticks(rotation=90)
# xlim_left=dates.date2num(datetime.datetime(2016,8,29,17,30,0))
# xlim_right=dates.date2num(datetime.datetime(2016,8,29,20,23,0))
# ax.set_xlim(xlim_left,xlim_right)
# ax.set_ylim(150000,350000)
# plt.grid()
# plt.title("Airglow's center height 160829")
# plt.savefig("exam_pumping_scheme_160829.png")
# # plt.close()
# plt.show()


# In[ ]:



