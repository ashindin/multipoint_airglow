
# coding: utf-8

# In[71]:

import time
import os, sys, inspect
import datetime
import numpy as np
import numexpr as ne

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import dates

import scipy.optimize as so


# In[2]:

from astropy.io import fits


# In[3]:

cmd_subfolder1 = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../astrometric_calibration")))
cmd_subfolder2 = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../modelling")))
if cmd_subfolder1 not in sys.path:
    sys.path.insert(0, cmd_subfolder1)
if cmd_subfolder2 not in sys.path:
    sys.path.insert(0, cmd_subfolder2)


# In[4]:

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)


# In[5]:

from arc_module import *
from tan_module import *


# In[6]:

#from model_sphere_fun import *
#from model_spheroid_fun import *
#from model_ellipsoid_fun import *
from model_spheroid_incl_fun import *
#from model_drop_fun import *

# In[86]:

# In[8]:

spath="./sphere/"
if not os.path.exists(spath):
    os.makedirs(spath)
spath="./spheroid/"
if not os.path.exists(spath):
    os.makedirs(spath)
spath="./ellipsoid/"
if not os.path.exists(spath):
    os.makedirs(spath)
spath="./spheroid_incl/"
if not os.path.exists(spath):
    os.makedirs(spath)
spath="./drop/"
if not os.path.exists(spath):
    os.makedirs(spath)

# In[74]:

def simpson(fun,a,b,n,args=()):
    h=(b-a)/n
    
    i=1
    x=a    
    S=fun(x,*args)
    
    x+=h    
    S+=4*fun(x,*args)
    i+=2
    
    while  i<n:
        x+=h        
        S+=2*fun(x,*args)
        x+=h        
        S+=4*fun(x,*args)
        i+=2        
    
    return (S+fun(b,*args))*h/3

def normalize_args(mod_pos, a):
    b=[0,0,0,0,0,0,0,0]
    b[0]=(mod_pos[0]-55.87*np.pi/180)/((56.4-55.87)*np.pi/180) # model lat
    b[1]=(mod_pos[1]-45.55*np.pi/180)/((46.55-45.55)*np.pi/180) # model lon
    b[2]=(mod_pos[2]-150000.)/(350000.-150000.) 
    
    b[3]=(a[0]-250.)/2250. # model concentration    
    b[4]=(a[1]-5000.)/(25000.) # model radius1    
    b[5]=(a[2]-5000.)/(25000.) # model radius2
    
    b[6]=a[3]/2/np.pi # model angle 1
    b[7]=a[4]/2/np.pi # model angle 2
    
    return b

def denormalize_args(b):
    mod_pos=[0,0,0]
    
    a=[0,0,0,0,0]
    
    mod_pos[0]=b[0]*((56.4-55.87)*np.pi/180)+55.87*np.pi/180 # model lat
    mod_pos[1]=b[1]*((46.55-45.55)*np.pi/180)+45.55*np.pi/180 # model lon
    mod_pos[2]=b[2]*(350000.-150000.)+150000. # model h
    
    a[0]=b[3]*2250. + 250. # model concentration   
    a[1]=b[4]*(25000.)+5000. # model radius1
    a[2]=b[5]*(25000.)+5000. # model radius2
    
    a[3]=b[6]*2*np.pi
    a[4]=b[7]*2*np.pi
    
    return mod_pos, a

def comparing_fun(b,model_fun,img1,img2,img1_ALT,img1_AZ,img2_ALT,img2_AZ,cam1_pos,cam2_pos):
    
    mod_pos, a = denormalize_args(b)
    
    m1=simpson(model_fun,200000.,400000.,100,(img1_ALT,img1_AZ, a, mod_pos,cam1_pos))
    m2=simpson(model_fun,150000.,350000.,100,(img2_ALT,img2_AZ, a ,mod_pos,cam2_pos))
    
#     m1=model_fun(a,r_s1c,lat_s1c,lon_s1c,az_pix_s1c,tet_pix_s1c)
#     m2=model_fun(a,r_keo,lat_keo,lon_keo,az_pix_keo,tet_pix_keo)
#     ret=np.sum((im1-m1)**2)+np.sum((im2-m2)**2)
    ret=0.12056327160493827*np.sum((img1-m1)**2)+(np.sum((img2-m2)**2))*0.25
    return ret

def represent_ans(x):
    return ' '.join(("{0:.3f}".format(x[0][0]*180/np.pi), "{0:.3f}".format(x[0][1]*180/np.pi), str(int(x[0][2])), str(int(x[1][0])), str(int(x[1][1]))))  


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

def plotting(couple_fn,img1,img2,m1,m2,img1_data_obs,img2_data_obs,img1_exptime,img2_exptime,out_path,out_fn):
    # plotting
    img1_data_obs_x=dates.date2num(img1_data_obs)
    img2_data_obs_x=dates.date2num(img2_data_obs)
    img1_data_end_x=dates.date2num(img1_data_obs+datetime.timedelta(seconds=img1_exptime))
    img2_data_end_x=dates.date2num(img2_data_obs+datetime.timedelta(seconds=img2_exptime))
    img1_spcalo=fits.getval(couple_fn,'SPCAL-O',0)
    img2_spcalo=fits.getval(couple_fn,'SPCAL-O',1)

    pumping_scheme_str_list=[]
    s_i=0
    while True:
        try:
            temp=fits.getval(couple_fn,'P-SCH-'+str(s_i),0)
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

    vmin=-15
    vmax=15
    fig=plt.figure(figsize=(12.8,7.2))
    fig.set_size_inches(12.8, 7.2)
    dx=0.05
    ax1=plt.axes(position=[0.035000000000000003+dx/2, 0.26, 0.4, 0.7])    
    pcm1=plt.pcolormesh(img1,vmin=vmin,vmax=vmax)    
    
    ret=0
    try:
        CS1 = plt.contour(m1, 2, colors='k')    
        plt.clabel(CS1, fontsize=9, inline=1,fmt='%1.1f')
    except:
        ret=1

    plt.title("CAM1: "+ img1_data_obs.strftime("%Y-%m-%dT%H:%M:%S.%f"), loc='left')
    plt.title("SP_OFFSET: " + "{:4.2f}".format(img1_spcalo),loc='right')
    ax1.set_ylim(101,0)
    ax1.set_xlim(101,0)
    plt.axis('equal')
    ax1_cb=plt.axes(position=[0.035000000000000003+dx/2+0.4+0.008, 0.26, 0.015, 0.7])
    plt.colorbar(pcm1,ax1_cb)

    ax2=plt.axes(position=[0.485+0.035000000000000003+dx/2, 0.26, 0.4, 0.7])
    pcm2=plt.pcolormesh(np.fliplr(img2),vmin=vmin,vmax=vmax)
    try:
        CS2 = plt.contour(np.fliplr(m2), 2, colors='k')
        plt.clabel(CS2, fontsize=9, inline=1,fmt='%1.1f')
    except:
        ret=1
    plt.title("CAM2: "+ img2_data_obs.strftime("%Y-%m-%dT%H:%M:%S.%f"),loc='left')
    plt.title("SP_OFFSET: " + "{:4.2f}".format(img2_spcalo),loc='right')    
    ax2.set_ylim(img2.shape[0],0)
#     ax2.set_xlim(0,img2.shape[1])
    ax2.set_xlim(img2.shape[1],0)
    plt.axis('equal')
    ax2_cb=plt.axes(position=[0.485+0.035000000000000003+dx/2+0.4+0.008, 0.26, 0.015, 0.7])
    plt.colorbar(pcm2,ax2_cb)

    ax3=plt.axes(position=[0.035000000000000003+dx/2, 0.05, 0.907, 0.17])
    xlim_left=img1_data_obs_x-12/60/24 #
    xlim_right=img1_data_obs_x+12/60/24 #
    ax3.set_xlim(xlim_left,xlim_right) 
    ax3.set_ylim(0,1.)
    ax3.xaxis.set_major_locator(dates.MinuteLocator(byminute=range(1,61,3)))
    ax3.xaxis.set_minor_locator(dates.MinuteLocator(byminute=range(4,64,3)))
    xfmt = dates.DateFormatter('%H:%M')
    ax3.xaxis.set_major_formatter(xfmt)
    for i in range(len(clean_pumping)):
        plt.plot(x_dates_clean[i],clean_pumping[i]*0.9,'r')

        plt.plot([img1_data_obs_x, img1_data_end_x],[0.95,0.95],'b',lw=3) #
        plt.plot([img2_data_obs_x, img2_data_end_x],[0.92,0.92],'g',lw=3) #
        plt.grid()
        plt.gca().yaxis.set_major_locator(plt.NullLocator())

#     plt.show()
    plt.savefig(out_path+out_fn[0:-4]+".png", dpi = 100)    
    plt.close()

    return ret
# In[82]:

def inv_problem_solve(couple_fn,model_fun,out_path):
    t=time.time()
    img1=fits.getdata(couple_fn,0)
    img2=fits.getdata(couple_fn,1)
    
    img1_data_obs=datetime.datetime.strptime(fits.getval(couple_fn,'DATE-OBS',0),"%Y-%m-%dT%H:%M:%S.%f")
    img2_data_obs=datetime.datetime.strptime(fits.getval(couple_fn,'DATE-OBS',1),"%Y-%m-%dT%H:%M:%S.%f")
    
    img1_exptime=fits.getval(couple_fn,'EXPTIME',0)
    img2_exptime=fits.getval(couple_fn,'EXPTIME',1)
    
    out_fn=(img1_data_obs+datetime.timedelta(seconds=img1_exptime/2)).strftime('%y%m%d_%H%M%S_%f_'+out_path.split('/')[1]+'.dat')
    
    prev_mod_type='spheroid'
    fn_previous='./'+prev_mod_type+'/'+(img1_data_obs+datetime.timedelta(seconds=img1_exptime/2)).strftime('%y%m%d_%H%M%S_%f_'+prev_mod_type+'.dat')
    fid=open(fn_previous,'r')
    fid.readline()
    success_str=fid.readline()
    status_str=fid.readline()
    fid.close()
    
    if success_str!="SUCCESS = True\n" or status_str!="STATUS = 0\n":
        return 1
    
    fid=open(fn_previous,'r')
    fid.readline()
    fid.readline()
    fid.readline()
    fid.readline()
    fid.readline()
    fid.readline()
    args_prev_str=fid.readline()[0:-1]
    args_prev=[float(ap) for ap in args_prev_str.split(' ')]
    fid.close()
    
    img1_lat=fits.getval(couple_fn,'LAT-DEG',0)*np.pi/180
    img2_lat=fits.getval(couple_fn,'LAT-DEG',1)*np.pi/180
    img1_lon=fits.getval(couple_fn,'LON-DEG',0)*np.pi/180
    img2_lon=fits.getval(couple_fn,'LON-DEG',1)*np.pi/180
    img1_hei=fits.getval(couple_fn,'HEI_M',0)
    img2_hei=fits.getval(couple_fn,'HEI_M',1)
    
    cam1_pos=(img1_lat,img1_lon,img1_hei)
    cam2_pos=(img2_lat,img2_lon,img2_hei)
    
    img1_az0=fits.getval(couple_fn,'AMCAL-A0',0)
    img2_az0=fits.getval(couple_fn,'AMCAL-A0',1)
    
    img1_proj=fits.getval(couple_fn,'PROJ-TYP',0)
    img2_proj=fits.getval(couple_fn,'PROJ-TYP',1)
    
    img1_alt0=fits.getval(couple_fn,'AMCAL-H0',0)
    img2_alt0=fits.getval(couple_fn,'AMCAL-H0',1)
    
    img1_a=fits.getval(couple_fn,'AMCAL-A',0)[1:-1].split()
    img1_a=[float(img1_a[i]) for i in range(len(img1_a))]
    img1_b=fits.getval(couple_fn,'AMCAL-B',0)[1:-1].split()
    img1_b=[float(img1_b[i]) for i in range(len(img1_b))]
    img1_c=fits.getval(couple_fn,'AMCAL-C',0)[1:-1].split()
    img1_c=[float(img1_c[i]) for i in range(len(img1_c))]
    img1_d=fits.getval(couple_fn,'AMCAL-D',0)[1:-1].split()
    img1_d=[float(img1_d[i]) for i in range(len(img1_d))]
    
    img2_a=fits.getval(couple_fn,'AMCAL-A',1)[1:-1].split()
    img2_a=[float(img2_a[i]) for i in range(len(img2_a))]
    img2_b=fits.getval(couple_fn,'AMCAL-B',1)[1:-1].split()
    img2_b=[float(img2_b[i]) for i in range(len(img2_b))]
    img2_c=fits.getval(couple_fn,'AMCAL-C',1)[1:-1].split()
    img2_c=[float(img2_c[i]) for i in range(len(img2_c))]
    img2_d=fits.getval(couple_fn,'AMCAL-D',1)[1:-1].split()
    img2_d=[float(img2_d[i]) for i in range(len(img2_d))]
    
    img1_YPIX, img1_XPIX = np.mgrid[1:img1.shape[0]+1, 1:img1.shape[1]+1]
    img2_YPIX, img2_XPIX = np.mgrid[1:img2.shape[0]+1, 1:img2.shape[1]+1]
    
    if img1_proj=='ARC':
        img1_AZ,img1_ALT=arc_pix2hor(img1_XPIX,img1_YPIX,img1_az0,img1_alt0,img1_a,img1_b)        
        img1_AZ = img1_AZ[325:426,175:276]
        img1_ALT = img1_ALT[325:426,175:276]
        img1=img1[325:426,175:276]
    elif img1_proj=='TAN':
        img1_AZ,img1_ALT=tan_pix2hor(img1_XPIX,img1_YPIX,img1_az0,img1_alt0,img1_a,img1_b)
        
    if img2_proj=='ARC':
        img2_AZ,img2_ALT=arc_pix2hor(img2_XPIX,img2_YPIX,img2_az0,img2_alt0,img2_a,img2_b)
        img2_AZ = img2_AZ[325:426,175:276]
        img2_ALT = img2_ALT[325:426,175:276]
        img2=img2[325:426,175:276]
    elif img2_proj=='TAN':
        img2_AZ,img2_ALT=tan_pix2hor(img2_XPIX,img2_YPIX,img2_az0,img2_alt0,img2_a,img2_b)    
    
    mod_pos= (args_prev[0]*np.pi/180, args_prev[1]*np.pi/180, args_prev[2])
    mod_args= (args_prev[3], args_prev[4], args_prev[5], 0., 0.)
    
    res = so.minimize(comparing_fun, normalize_args(mod_pos,mod_args), (model_fun,img1,img2,img1_ALT,img1_AZ,img2_ALT,img2_AZ,cam1_pos,cam2_pos),                       method='Nelder-Mead',options=
                      {'return_all':True,'maxiter':4000, 'maxfev':4000,'xatol': 0.0001,'fatol': 0.0001})
    
    if res.status==0 and res.success==True:
        fid=open(out_path+out_fn,'w')

        
        denorm_res=denormalize_args(res.x)
        denorm_x=denormalize_args(res.x)
        denorm_x[0][0]*=180/np.pi
        denorm_x[0][1]*=180/np.pi        
        denorm_x_str_list1=[str(dx) for dx in denorm_x[0]]
        denorm_x_str_list2=[str(dx) for dx in denorm_x[1]]
        denorm_x_str_list=denorm_x_str_list1+denorm_x_str_list2        
        denorm_x_str=' '.join(denorm_x_str_list)
        
        fid.write('# '+couple_fn.split('/')[-1]+'\n')
        
        if all([dx>0 for dx in denorm_x[0]])>0 and all([dx>0 for dx in denorm_x[1][0:-2]])>0: 

            m1=simpson(model_fun,200000.,400000.,100,(img1_ALT,img1_AZ, denorm_res[1], denorm_res[0],cam1_pos))
            m2=simpson(model_fun,150000.,350000.,100,(img2_ALT,img2_AZ, denorm_res[1], denorm_res[0],cam2_pos))
            plot_ret = plotting(couple_fn,img1,img2,m1,m2,img1_data_obs,img2_data_obs,img1_exptime,img2_exptime,out_path,out_fn)
            if plot_ret==0:                
                fid.write('SUCCESS = '+str(res.success)+'\n')
            else:
                fid.write('SUCCESS = '+'False'+'\n')

        else:
            fid.write('SUCCESS = '+'False'+'\n')         
            
        fid.write('STATUS = '+str(res.status)+'\n')
        fid.write('NFEV = '+str(res.nfev)+'\n')
        fid.write('NIT = '+str(res.nit)+'\n')
        fid.write('FUN = '+str(res.fun)+'\n')
                
        fid.write(denorm_x_str+'\n\n')
        
        for i in range(len(res.allvecs)):  
            denorm_x=denormalize_args(res.allvecs[i])
            denorm_x[0][0]*=180/np.pi
            denorm_x[0][1]*=180/np.pi        
            denorm_x_str_list1=[str(dx) for dx in denorm_x[0]]
            denorm_x_str_list2=[str(dx) for dx in denorm_x[1]]
            denorm_x_str_list=denorm_x_str_list1+denorm_x_str_list2
            denorm_x_str=' '.join(denorm_x_str_list)
            fid.write(denorm_x_str+'\n')         
        fid.close()
        
    else:
        fid=open(out_path+out_fn,'w')
        fid.write('# '+couple_fn.split('/')[-1]+'\n')        
        fid.write('SUCCESS = '+str(res.success)+'\n')        
        fid.write('STATUS = '+str(res.status)+'\n')
        fid.write('NFEV = '+str(res.nfev)+'\n')
        fid.write('NIT = '+str(res.nit)+'\n\n')
        
        for i in range(len(res.allvecs)):  
            denorm_x=denormalize_args(res.allvecs[i])
            denorm_x[0][0]*=180/np.pi
            denorm_x[0][1]*=180/np.pi        
            denorm_x_str_list1=[str(dx) for dx in denorm_x[0]]
            denorm_x_str_list2=[str(dx) for dx in denorm_x[1]]
            denorm_x_str_list=denorm_x_str_list1+denorm_x_str_list2
            denorm_x_str=' '.join(denorm_x_str_list)
            fid.write(denorm_x_str+'\n')   
        
        fid.close()
    
    elapsed=time.time()-t
    print(couple_fn.split('/')[-1],'-','SUCCESS = '+str(res.success),'STATUS = '+str(res.status),
          'NFEV = '+str(res.nfev),'NIT = '+str(res.nit),'-', elapsed,'sec')
    return 0
# t=time.time()
# res = inv_problem_solve(couple_fn,sphere_fun,((56.1434444*np.pi/180,46.0991056*np.pi/180,230000.),(850,19000.)),'./sphere/')
# elapsed=time.time()-t
# print(elapsed)
# res


# In[ ]:

def inv_problem_solve_spath(spath):
    
    couples_fit_filenames=[]
    for fn in next(os.walk(spath))[2]:
        if fn[-4:]=='fits':
            couples_fit_filenames.append(spath+fn)
    couples_fit_filenames=sorted(couples_fit_filenames)
    
    #for i in range(len(couples_fit_filenames)):
    for i in range(0,1):
        #sys.stdout.write('\r')
        #sys.stdout.write("Processing frame "+str(i+1)+"/"+str(len(couples_fit_filenames)))
        #sys.stdout.flush()

        couple_fn=couples_fit_filenames[i]
        inv_problem_solve(couple_fn,spheroid_incl_fun,'./spheroid_incl/')
    return None


# In[ ]:

inv_problem_solve_spath("../image_processing/couples140824/")
#inv_problem_solve_spath("../image_processing/couples140826/")
#inv_problem_solve_spath("../image_processing/couples160829/")

