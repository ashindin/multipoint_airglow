
# coding: utf-8

# In[1]:

import datetime
import os, sys, inspect
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import dates


# In[2]:

cmd_subfolder1 = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../astrometric_calibration")))
cmd_subfolder2 = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../modelling")))
if cmd_subfolder1 not in sys.path:
    sys.path.insert(0, cmd_subfolder1)
if cmd_subfolder2 not in sys.path:
    sys.path.insert(0, cmd_subfolder2)


# In[3]:

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from arc_module import *
from tan_module import *


# In[4]:

from model_drop_fun import *


# In[5]:

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


# In[6]:

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


# In[7]:

def normalize_args(mod_pos, a):
    b=[0,0,0,0,0,0,0,0]
    b[0]=(mod_pos[0]-55.87*np.pi/180)/((56.4-55.87)*np.pi/180) # model lat
    b[1]=(mod_pos[1]-45.55*np.pi/180)/((46.55-45.55)*np.pi/180) # model lon
    b[2]=(mod_pos[2]-150000.)/(350000.-150000.) 
    
    b[3]=(a[0]-250.)/2250. # model concentration    
    b[4]=(a[1]-5000.)/(25000.) # model radius1    
    b[5]=(a[2]-5000.)/(25000.) # model radius2
    
    b[6]=(a[3]-5000.)/(25000.) # model height
    b[7]=a[4] # drop parameter
    
    return b

def denormalize_args(b):
    mod_pos=[0,0,0]
    
    a=[0,0,0,0,0]
    
    mod_pos[0]=b[0]*((56.4-55.87)*np.pi/180)+55.87*np.pi/180 # model lat
    mod_pos[1]=b[1]*((46.55-45.55)*np.pi/180)+45.55*np.pi/180 # model lon
    mod_pos[2]=b[2]*(350000.-150000.)+150000. # model h
    
    a[0]=b[3]*2250. + 250. # model concentration   
    a[1]=b[4]*(25000.)+5000. # model radius1
    a[2]=b[5]*(25000.)+5000. # model radius1
    a[3]=b[6]*(25000.)+5000. # model height
    a[4]=b[7] # drop parameter
    
    return mod_pos, a


# In[8]:

couple_fn='../image_processing/couples140824/couple_140824_179.fits'
dat_fn='../solving/drop/140824_193947_000000_drop.dat'

# couple_fn='../image_processing/couples140824/couple_140824_154.fits'
# dat_fn='../solving/drop/140824_192551_000000_drop.dat'

keo_solve_pars_fn="../astrometric_calibration/keo_140824_solve_manual.pars"
s1c_solve_pars_fn="../astrometric_calibration/s1c_140824_solve.pars"


# In[9]:

img1_err,img1_az0,img1_alt0,img1_a,img1_b,img1_c,img1_d=get_solve_pars(keo_solve_pars_fn)
img2_err,img2_az0,img2_alt0,img2_a,img2_b,img2_c,img2_d=get_solve_pars(s1c_solve_pars_fn)


# In[10]:

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


# In[11]:

img1_lat=fits.getval(couple_fn,'LAT-DEG',0)*np.pi/180
img2_lat=fits.getval(couple_fn,'LAT-DEG',1)*np.pi/180
img1_lon=fits.getval(couple_fn,'LON-DEG',0)*np.pi/180
img2_lon=fits.getval(couple_fn,'LON-DEG',1)*np.pi/180
img1_hei=fits.getval(couple_fn,'HEI_M',0)
img2_hei=fits.getval(couple_fn,'HEI_M',1)

cam1_pos=(img1_lat,img1_lon,img1_hei)
cam2_pos=(img2_lat,img2_lon,img2_hei)


# In[12]:

cam1_pos


# In[13]:

img1=fits.getdata(couple_fn,0)
# img1 = img1[325:426,175:276]
img2=fits.getdata(couple_fn,1)

img1_YPIX, img1_XPIX = np.mgrid[1:img1.shape[0]+1, 1:img1.shape[1]+1]
img2_YPIX, img2_XPIX = np.mgrid[1:img2.shape[0]+1, 1:img2.shape[1]+1]

img1_proj=fits.getval(couple_fn,'PROJ-TYP',0)
img2_proj=fits.getval(couple_fn,'PROJ-TYP',1)
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


# In[14]:

fun, mod_pos, pars = get_model_pars(dat_fn)
# print(np.array(mod_pos[0:2])*180/np.pi,mod_pos[2])
# print(pars)


# In[ ]:

m1=simpson(drop_fun,200000.,400000.,1000,(img1_ALT,img1_AZ, pars, mod_pos,cam1_pos))
m2=simpson(drop_fun,150000.,350000.,1000,(img2_ALT,img2_AZ, pars, mod_pos,cam2_pos))

vmin=0
vmax=15

fig=plt.figure(figsize=(6.69,3.4))

ax1=plt.axes(position=[0.055, 0.1, 0.4, 0.79])   
ax2=plt.axes(position=[0.55-0.035, 0.1, 0.4, 0.79]) 

plt.sca(ax1)
pcm1=plt.pcolormesh(img1,vmin=vmin,vmax=vmax,cmap='jet') 
# pcm1=plt.pcolormesh(m1,vmin=vmin,vmax=vmax,cmap='jet') 
CS1 = plt.contour(m1, 2, colors='k',linewidths=3)    
plt.clabel(CS1, fontsize=10, inline=1,fmt='%1.1f')
pcm1.set_edgecolor('face')
plt.ylim(100,0)
plt.xlim(100,0)
plt.axis('equal')

plt.sca(ax2)
pcm2=plt.pcolormesh(img2,vmin=vmin,vmax=vmax,cmap='jet')
# pcm2=plt.pcolormesh(m2,vmin=vmin,vmax=vmax,cmap='jet')
CS2 = plt.contour(m2, 2, colors='k',linewidths=3)    
plt.clabel(CS2, fontsize=10, inline=1,fmt='%1.1f')
pcm2.set_edgecolor('face')
plt.ylim(288,0)
plt.xlim(0,288)
plt.axis('equal')

ax2_cb=plt.axes(position=[0.55+0.43+0.008-0.065, 0.1, 0.015, 0.79])
plt.colorbar(pcm2,ax2_cb)  
plt.savefig('fig7.png',dpi=300)
plt.savefig('fig7.pdf',dpi=300)
plt.savefig('fig7.eps',dpi=300)
plt.show()

