
# coding: utf-8

# In[3]:

import os, sys, inspect
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.integrate as si
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../modelling")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
from model_drop_fun import *
import time 
import numexpr as ne
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../astrometric_calibration")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
from tan_module import *


# In[4]:

from sympy import *


# In[5]:

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
    return err, az0,alt0,a,b,c,d

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
# simpson(sphere_fun,0.,500000.,120,(E,A,aa,mod_pos,cam_pos))


# In[6]:

mod_pos=(np.pi/2,0.,250000.)
cam_pos=(np.pi/2,0.,0.)
E=np.pi/2
A=0.
R=np.arange(0,500000.,1000.)


# In[7]:

solve_pars_fname="../astrometric_calibration/s1c_140824_solve.pars"
err,az0,alt0,a,b,c,d=get_solve_pars(solve_pars_fname)
YPIX, XPIX = np.mgrid[1:288+1, 1:288+1]
AZ,ALT=tan_pix2hor(XPIX,YPIX,az0,alt0,a,b)


# In[23]:

z_axe=np.linspace(0.,500000.,5000)
shift_axe=np.linspace(0,1,5)

# print(shift_axe)
fig=plt.figure(figsize=(6,3))
ax=plt.axes(position=[0.1,0.17,0.85,0.8])
lines=[]
leg=['$p_5=2, p_6=6$','$p_5=4, p_6=6$','$p_5=6, p_6=6$','$p_5=6, p_6=4$','$p_5=6, p_6=2$']
for shift in shift_axe:
    i=0
    b_axe=drop_fun(z_axe,np.pi/2,0.,(1000.,10000.*2**.5, 10000.*2**.5, 70000., shift),mod_pos,cam_pos)/(5.627/10**7)
    line, = plt.plot(z_axe-250000,b_axe)
    lines.append(line)
    i+=1
plt.grid()
plt.xlim(250000-250000-5000,320000-250000+35000)
plt.ylim(0,1100)
plt.yticks([0,500,1000], ['0','','$p_1$'])
plt.xticks([0,35000,70000], ['0','','$p_4$'])
plt.xlabel('$z$')
plt.ylabel('$n(z)$')
plt.legend(lines,leg)
# plt.show()
plt.savefig("fig6.png",dpi=300)
plt.savefig("fig6.eps",dpi=300)
plt.savefig("fig6.pdf",dpi=300)
# plt.show()

