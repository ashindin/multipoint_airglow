import datetime
import time
import subprocess
import os
import platform
from pathlib import Path
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import scipy.optimize as so
import numpy as np
import matplotlib.pyplot as plt
os_name=platform.system()

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

from tan_module import *
from save_solve_pars import *

win_com_prefix='bash --login -c "('
win_com_postfix=')"'
s1c_solve_pars='--scale-units arcsecperpix --scale-low 250 --scale-high 270'

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

def image2xy(fname, out_dir):
    os_name=platform.system()
#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)
    name=fname.split('/')[-1]
    out_name=out_dir + '/' + name.split('.')[-2] + '.axy'
    com_line='/usr/local/astrometry/bin/image2xy -O -o ' + out_name +' ' + fname
    if os_name=='Windows':
        com_line=win_com_prefix+com_line+win_com_postfix
#     print(com_line)
    return  os.system(com_line), out_name
# image2xy(fname, spath)

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
    
    com_line='cd ' + spath  + ' && /usr/local/astrometry/bin/solve-field ' + axy_fname.split('/')[-1] +' --continue -D ' + spath + ' ' + solve_pars + ' --cpulimit 2 --no-plots -M none -S none -B none -W none'
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
solve_fname="solve_field_altaz_140824_s1c.dat"
save_fname="s1c_140824_solve.pars"
save_solve_data(fit_path,solve_fname)
save_med_solve_pars(save_fname)

fit_path="../data/140826/s1c"
solve_fname="solve_field_altaz_140826_s1c.dat"
save_fname="s1c_140826_solve.pars"
save_solve_data(fit_path,solve_fname)
save_med_solve_pars(save_fname)
