
# coding: utf-8

# In[1]:

import datetime
import time
import subprocess
import os
import platform
from pathlib import Path
# import fitsio
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import scipy.optimize as so
import numpy as np
import matplotlib.pyplot as plt
os_name=platform.system()


# In[2]:

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)


# In[3]:

win_com_prefix='bash --login -c "('
win_com_postfix=')"'


# In[4]:

s1c_solve_pars='--scale-units arcsecperpix --scale-low 250 --scale-high 270'
keo_solve_pars=''


# In[5]:

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


# In[6]:

def keo_get_date_obs(filename):
    fid_fit=fits.open(filename);
    date_obs_str=fid_fit[0].header["DATE-OBS"]  
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S")


# In[7]:

def tan_calc_pix2st_coefs(arg, AZ, ALT, x, y):
    
    az0=arg[0]
    alt0=arg[1]
    
    on=np.ones(len(x))
    
    X=-(180/np.pi)*np.cos(ALT)*np.sin(AZ-az0)/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    Y=(180/np.pi)*(np.sin(alt0)*np.cos(ALT)*np.cos(AZ-az0)-np.cos(alt0)*np.sin(ALT))/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    
    matA=np.vstack((on,X,Y)).T
    matB=np.vstack((on,x,y)).T
    
    imatBBB=np.dot(np.linalg.inv(np.dot(matB.T,matB)),matB.T)
    imatAAA=np.dot(np.linalg.inv(np.dot(matA.T,matA)),matA.T)
    
    a=np.dot(imatBBB,X).T
    b=np.dot(imatBBB,Y).T
    c=np.dot(imatAAA,x).T
    d=np.dot(imatAAA,y).T
    
    #x_new=c[0]+c[1]*X+c[2]*Y;
    #y_new=d[0]+d[1]*X+d[2]*Y;
    
#     sigma=np.sqrt(np.mean((x-x_new)**2+(y-y_new)**2));
    
    return a,b,c,d


# In[8]:

def arc_calc_pix2st_coefs(arg, AZ, ALT, x, y):
    
    az0=arg[0]
    alt0=arg[1]
    
    on=np.ones(len(x))
    
    ALT_I=np.arcsin(np.sin(ALT)*np.sin(alt0)+np.cos(ALT)*np.cos(alt0)*np.cos(AZ-az0))
    AZ_I=np.arctan2(-np.cos(ALT)*np.sin(AZ-az0),np.sin(ALT)*np.cos(alt0)-np.cos(ALT)*np.sin(alt0)*np.cos(AZ-az0))
    
    X=(90-ALT_I*180/np.pi)*np.sin(AZ_I)    
    Y=(90-ALT_I*180/np.pi)*np.cos(AZ_I)
    
    matA=np.vstack((on,X,Y)).T
    matB=np.vstack((on,x,y)).T
    
    imatBBB=np.dot(np.linalg.inv(np.dot(matB.T,matB)),matB.T)
    imatAAA=np.dot(np.linalg.inv(np.dot(matA.T,matA)),matA.T)
    
    a=np.dot(imatBBB,X).T
    b=np.dot(imatBBB,Y).T
    c=np.dot(imatAAA,x).T
    d=np.dot(imatAAA,y).T
    
    #x_new=c[0]+c[1]*X+c[2]*Y;
    #y_new=d[0]+d[1]*X+d[2]*Y;
    
#     sigma=np.sqrt(np.mean((x-x_new)**2+(y-y_new)**2));
    
    return a,b,c,d


# In[9]:

def tan_calc_pix2st_coefs_discrep(arg, AZ, ALT, x, y):
    
    az0=arg[0]
    alt0=arg[1]
    
    on=np.ones(len(x))
    
    X=-(180/np.pi)*np.cos(ALT)*np.sin(AZ-az0)/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    Y=(180/np.pi)*(np.sin(alt0)*np.cos(ALT)*np.cos(AZ-az0)-np.cos(alt0)*np.sin(ALT))/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    
    matA=np.vstack((on,X,Y)).T
    matB=np.vstack((on,x,y)).T
    
    imatBBB=np.dot(np.linalg.inv(np.dot(matB.T,matB)),matB.T)
    imatAAA=np.dot(np.linalg.inv(np.dot(matA.T,matA)),matA.T)
    
    a=np.dot(imatBBB,X).T
    b=np.dot(imatBBB,Y).T
    c=np.dot(imatAAA,x).T
    d=np.dot(imatAAA,y).T
    
    x_new=c[0]+c[1]*X+c[2]*Y;
    y_new=d[0]+d[1]*X+d[2]*Y;
    
    sigma=np.sqrt(np.mean((x-x_new)**2+(y-y_new)**2));
    
    return sigma


# In[10]:

def arc_calc_pix2st_coefs_discrep(arg, AZ, ALT, x, y):
    
    az0=arg[0]
    alt0=arg[1]
    
    on=np.ones(len(x))
    
    ALT_I=np.arcsin(np.sin(ALT)*np.sin(alt0)+np.cos(ALT)*np.cos(alt0)*np.cos(AZ-az0))
    AZ_I=np.arctan2(-np.cos(ALT)*np.sin(AZ-az0),np.sin(ALT)*np.cos(alt0)-np.cos(ALT)*np.sin(alt0)*np.cos(AZ-az0))
    
    X=(90-ALT_I*180/np.pi)*np.sin(AZ_I)    
    Y=(90-ALT_I*180/np.pi)*np.cos(AZ_I)
        
    matA=np.vstack((on,X,Y)).T
    matB=np.vstack((on,x,y)).T
    
    imatBBB=np.dot(np.linalg.inv(np.dot(matB.T,matB)),matB.T)
    imatAAA=np.dot(np.linalg.inv(np.dot(matA.T,matA)),matA.T)
    
    a=np.dot(imatBBB,X).T
    b=np.dot(imatBBB,Y).T
    c=np.dot(imatAAA,x).T
    d=np.dot(imatAAA,y).T
    
    x_new=c[0]+c[1]*X+c[2]*Y;
    y_new=d[0]+d[1]*X+d[2]*Y;
    
    sigma=np.sqrt(np.mean((x-x_new)**2+(y-y_new)**2));
    
    return sigma


# In[11]:

def tan_hor2pix(AZ,ALT,az0,alt0,c,d):
    X=-(180/np.pi)*np.cos(ALT)*np.sin(AZ-az0)/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    Y=(180/np.pi)*(np.sin(alt0)*np.cos(ALT)*np.cos(AZ-az0)-np.cos(alt0)*np.sin(ALT))/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    X_PIX=c[0]+c[1]*X+c[2]*Y;
    Y_PIX=d[0]+d[1]*X+d[2]*Y;
    return X_PIX, Y_PIX;


# In[12]:

def arc_hor2pix(AZ,ALT,az0,alt0,c,d):
    ALT_I=np.arcsin(np.sin(ALT)*np.sin(alt0)+np.cos(ALT)*np.cos(alt0)*np.cos(AZ-az0))
    AZ_I=np.arctan2(-np.cos(ALT)*np.sin(AZ-az0),np.sin(ALT)*np.cos(alt0)-np.cos(ALT)*np.sin(alt0)*np.cos(AZ-az0))
    X=(90-ALT_I*180/np.pi)*np.sin(AZ_I)    
    Y=(90-ALT_I*180/np.pi)*np.cos(AZ_I)
    
    X_PIX=c[0]+c[1]*X+c[2]*Y;
    Y_PIX=d[0]+d[1]*X+d[2]*Y;
    return X_PIX, Y_PIX;


# In[13]:

def tan_pix2hor(x,y,az0,alt0,a,b):
    X=a[0]+a[1]*x+a[2]*y
    Y=b[0]+b[1]*x+b[2]*y
    
    AZ_I=np.arctan2(X,-Y)
    R=np.sqrt(X**2+Y**2)    
    
    ALT_I=np.arctan(180/np.pi/R) # TAN PROJECTION
    
    AZ=az0+np.arctan2(-np.cos(ALT_I)*np.sin(AZ_I),
    np.sin(ALT_I)*np.cos(alt0)-np.cos(ALT_I)*np.sin(alt0)*np.cos(AZ_I))
    
    ALT=np.arcsin(np.sin(ALT_I)*np.sin(alt0)+np.cos(ALT_I)*np.cos(alt0)*np.cos(AZ_I))    
    
    return AZ, ALT


# In[14]:

def arc_pix2hor(x,y,az0,alt0,a,b):
    X=a[0]+a[1]*x+a[2]*y
    Y=b[0]+b[1]*x+b[2]*y
    
    AZ_I=np.arctan2(X,-Y)
    R=np.sqrt(X**2+Y**2)
    
    ALT_I=(90-R)*np.pi/180 # ARC PROJECTION
    
    AZ=az0+np.arctan2(-np.cos(ALT_I)*np.sin(AZ_I),
    np.sin(ALT_I)*np.cos(alt0)-np.cos(ALT_I)*np.sin(alt0)*np.cos(AZ_I))
    
    ALT=np.arcsin(np.sin(ALT_I)*np.sin(alt0)+np.cos(ALT_I)*np.cos(alt0)*np.cos(AZ_I))    
    
    return AZ, ALT


# In[15]:

def image2xy(fname, out_dir):
    os_name=platform.system()
#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)
    name=fname.split('/')[-1]
    out_name=out_dir + '/' + name.split('.')[-2] + '.axy'
    com_line='image2xy -O -o ' + out_name +' ' + fname
    if os_name=='Windows':
        com_line=win_com_prefix+com_line+win_com_postfix
#     print(com_line)
    return  os.system(com_line), out_name
# image2xy(fname, spath)


# In[19]:

def s1c_solve_field_altaz(fname, solve_pars, get_date_obs_fun=s1c_get_date_obs,lat_deg=56.1501667,lon_deg=46.1050833,hei_m=183.):
    os_name=platform.system()
    
    fname=os.path.abspath(fname)
#     print(fname)
    path, file = os.path.split(fname)
    spath=path+"/.temp"
    fname_axy_abs=spath+"/"+file[:-4]+".axy"
    if not os.path.exists(spath):
        os.makedirs(spath)
        
    if os_name=='Windows':
        fname="/cygdrive/"+fname.replace(":","").replace("\\","/")
#     print(fname)
    path, file = os.path.split(fname)
    spath=path+"/.temp"
#     print(spath)    
    
    err_code, axy_fname = image2xy(fname, spath)
    
    if err_code!=0:
        return 1
    
    com_line='cd ' + spath  + ' && solve-field ' + axy_fname.split('/')[-1] +' --continue -D ' + spath + ' ' + solve_pars + ' --cpulimit 2 --no-plots -M none -S none -B none -W none'
    if os_name=='Windows':
        com_line=win_com_prefix+com_line+win_com_postfix
#     print(com_line)
    err_code=os.system(com_line)
    if err_code!=0:
        return 2
    
    site=EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg, height=hei_m*u.m)
    

    
    fname_xy=axy_fname[0:-4]+'-indx.xyls'
    fname_rd=axy_fname[0:-4]+'.rdls'
    
    if os_name=='Windows':
        fname=fname[10::]
        fname=fname[0]+':'+fname[1::]
        fname_xy=fname_xy[10::]
        fname_xy=fname_xy[0]+':'+fname_xy[1::]    
        fname_rd=fname_rd[10::]
        fname_rd=fname_rd[0]+':'+fname_rd[1::]
#     print(fname_xy)
#     print(fname_rd)
#     print(axy_fname)
    if Path(fname_rd).is_file()==False:
        return np.nan,np.nan,np.array([np.nan, np.nan, np.nan]),np.array([np.nan, np.nan, np.nan]), np.array([np.nan, np.nan, np.nan]),np.array([np.nan, np.nan, np.nan])
        

    date_obs=get_date_obs_fun(fname,ut_shift=-4)

    hdulist = fits.open(fname_xy)


    tbdata = hdulist[1].data

    X=tbdata.field('X')
    Y=tbdata.field('Y')

    hdulist.close()


    hdulist = fits.open(fname_rd)

    tbdata = hdulist[1].data

    RA=tbdata.field('RA')
    DEC=tbdata.field('DEC')
    hdulist.close()

    SC=SkyCoord(RA, DEC, frame='icrs', unit='deg');

    AA=SC.transform_to(AltAz(obstime=date_obs, location=site,temperature=0*u.deg_C,pressure=1013*u.hPa,
                                           relative_humidity=0.5,obswl=630.0*u.nm))

    AZ=AA.az.rad
    ALT=AA.alt.rad

    res=so.minimize(tan_calc_pix2st_coefs_discrep,(0,np.pi/2),(AZ,ALT,X,Y),method='Nelder-Mead')

    res_x=res.x
    if res_x[1]>np.pi/2:
        res_x[1]=np.pi-res_x[1]
        res_x[0]=res_x[0]+np.pi
    if res_x[0]<0:
        res_x[0]=res_x[0]+2*np.pi
    if res_x[0]>2*np.pi:
        res_x[0]=res_x[0]%(2*np.pi)

    az0=res_x[0];
    alt0=res_x[1];
    
    a,b,c,d=tan_calc_pix2st_coefs((az0,alt0),AZ,ALT,X,Y)
#     print(az0*180/np.pi,alt0*180/np.pi)
    az_c, alt_c = tan_pix2hor(144,144,az0,alt0,a,b)
    print("Central pixel direction (AZ, ALT in deg):",az_c*180/np.pi,alt_c*180/np.pi)       
    return az0,alt0,az_c,alt_c,a,b,c,d


fit_path="../data/140824/s1c"
fit_filenames=[fit_path+'/'+fn for fn in next(os.walk(fit_path))[2]]
res_fid=open('solve_field_altaz_140824.dat','w')
res_fid.write("# filename.fit az_c alt_c az0 alt0 a[0] a[1] a[2] b[0] b[1] b[2] c[0] c[1] c[2] d[0] d[1] d[2]\n")
for i in range(len(fit_filenames)):
    ret = s1c_solve_field_altaz(fit_filenames[i],s1c_solve_pars)
    if len(ret)==8:
        az0=ret[0]
        alt0=ret[1]
        az_c=ret[2]
        alt_c=ret[3]
        a=ret[4]
        b=ret[5]
        c=ret[6]
        d=ret[7]
        str_to_file=fit_filenames[i].split("/")[-1]+ " " + str(az_c) + " " +str(alt_c)+ " " +str(az0)+ " " +str(alt0)+ " "
        str_to_file+=" " + str(a[0]) + " " + str(a[1]) + " " + str(a[2])
        str_to_file+=" " + str(b[0]) + " " + str(b[1]) + " " + str(b[2])
        str_to_file+=" " + str(c[0]) + " " + str(c[1]) + " " + str(c[2])
        str_to_file+=" " + str(d[0]) + " " + str(d[1]) + " " + str(d[2]) + "\n"
        res_fid.write(str_to_file)
    else:
        print(ret)
        print(str(i)," ",fit_filenames[i].split("/")[-1]," ERROR")
res_fid.close()

fit_path="../data/140826/s1c"
fit_filenames=[fit_path+'/'+fn for fn in next(os.walk(fit_path))[2]]
res_fid=open('solve_field_altaz_140826.dat','w')
res_fid.write("# filename.fit az_c alt_c az0 alt0 a[0] a[1] a[2] b[0] b[1] b[2] c[0] c[1] c[2] d[0] d[1] d[2]\n")
for i in range(len(fit_filenames)):
    ret = s1c_solve_field_altaz(fit_filenames[i],s1c_solve_pars)
    if len(ret)==8:
        az0=ret[0]
        alt0=ret[1]
        az_c=ret[2]
        alt_c=ret[3]
        a=ret[4]
        b=ret[5]
        c=ret[6]
        d=ret[7]
        str_to_file=fit_filenames[i].split("/")[-1]+ " " + str(az_c) + " " +str(alt_c)+ " " +str(az0)+ " " +str(alt0)+ " "
        str_to_file+=" " + str(a[0]) + " " + str(a[1]) + " " + str(a[2])
        str_to_file+=" " + str(b[0]) + " " + str(b[1]) + " " + str(b[2])
        str_to_file+=" " + str(c[0]) + " " + str(c[1]) + " " + str(c[2])
        str_to_file+=" " + str(d[0]) + " " + str(d[1]) + " " + str(d[2]) + "\n"
        res_fid.write(str_to_file)
    else:
        print(ret)
        print(str(i)," ",fit_filenames[i].split("/")[-1]," ERROR")
res_fid.close()