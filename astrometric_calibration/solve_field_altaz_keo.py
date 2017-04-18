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

from arc_module import *
from save_solve_pars import *

win_com_prefix='bash --login -c "('
win_com_postfix=')"'

keo_solve_pars='--depth 11-20'

def keo_get_date_obs(filename):
    fid_fit=fits.open(filename);
    date_obs_str=fid_fit[0].header["DATE-OBS"]  
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S")
    
def keo_crop(fname,fname2):
    hdulist = fits.open(fname,ignore_missing_end=True)
    img=hdulist[0].data
    img2=img[255-44:255+45,255-44:255+45]
    hdulist[0].data=img2#-dark
    hdulist[0].header['NAXIS1']=89
    hdulist[0].header['NAXIS2']=89
    hdulist.writeto(fname2,overwrite=True)
    hdulist.close()
    return None

def keo_solve_field_altaz(fname, solve_pars, get_date_obs_fun=keo_get_date_obs,lat_deg=55.9305361,lon_deg=48.7444861,hei_m=91.):
    os_name=platform.system()
    
    fname=os.path.abspath(fname)
        
#     print(fname)
    path, file = os.path.split(fname)
    spath=path+"/.temp"
    if not os.path.exists(spath):
        os.makedirs(spath)
        
    crop_fname=spath+'/'+file.split('.')[-2]+'_crop.'+file.split('.')[-1]
    #print(crop_fname)
    keo_crop(fname,crop_fname)
    
    if not os.path.exists(spath):
        os.makedirs(spath)
        
    if os_name=='Windows':
        fname_cyg="/cygdrive/"+fname.replace(":","").replace("\\","/")
#     print(fname)
    path, file = os.path.split(fname)
    spath=path+"/.temp"
    
#     correct_keo_xy(axy_fname);
    crop_fname_cyg="/cygdrive/"+crop_fname.replace(":","").replace("\\","/")
#     com_line='cd ' + spath  + ' && solve-field ' + axy_fname.split('/')[-1] +' --continue -D ' + spath + ' ' + solve_pars + ' --cpulimit 2 --no-plots -M none -S none -B none -W none'
    com_line='/usr/local/astrometry/bin/solve-field ' + '--overwrite '  + solve_pars + ' --sigma 50 --crpix-center --cpulimit 2 --no-plots -M none -S none -B none -W none ' + crop_fname_cyg

    if os_name=='Windows':
        com_line=win_com_prefix+com_line+win_com_postfix
#     print(com_line)
    err_code=os.system(com_line)
    if err_code!=0:
        return 2
    
    site=EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg, height=hei_m*u.m)
        
    fname_xy=(crop_fname[0:-4]+'-indx.xyls').replace("\\","/")
    fname_rd=(crop_fname[0:-4]+'.rdls').replace("\\","/")
    
#     print(fname_xy)
#     print(fname_rd)    
    
    if Path(fname_rd).is_file()==False:
        return np.nan,np.nan,np.array([np.nan, np.nan, np.nan]),np.array([np.nan, np.nan, np.nan]), np.array([np.nan, np.nan, np.nan]),np.array([np.nan, np.nan, np.nan])
        
#     correct_keo_xyls(fname_xy)
        
    date_obs=get_date_obs_fun(fname)

    hdulist = fits.open(fname_xy)


    tbdata = hdulist[1].data

    X=tbdata.field('X')+211
    Y=tbdata.field('Y')+211

    hdulist.close()


    hdulist = fits.open(fname_rd)

    tbdata = hdulist[1].data

    RA=tbdata.field('RA')
    DEC=tbdata.field('DEC')


    SC=SkyCoord(RA, DEC, frame='icrs', unit='deg');

    AA=SC.transform_to(AltAz(obstime=date_obs, location=site,temperature=0*u.deg_C,pressure=1013*u.hPa,
                                           relative_humidity=0.5,obswl=630.0*u.nm))

    AZ=AA.az.rad
    ALT=AA.alt.rad

    res=so.minimize(arc_calc_pix2st_coefs_discrep,(0,np.pi/2),(AZ,ALT,X,Y),method='Nelder-Mead')

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
    
    a,b,c,d=arc_calc_pix2st_coefs((az0,alt0),AZ,ALT,X,Y)
#     print(az0*180/np.pi,alt0*180/np.pi)
    az_c, alt_c = arc_pix2hor(255,255,az0,alt0,a,b)
    print("Central pixel direction (AZ, ALT in deg):",az_c*180/np.pi,alt_c*180/np.pi)       
    return az0,alt0,az_c,alt_c,a,b,c,d

def save_solve_data(fit_path,solve_fname):
    fit_filenames=[fit_path+'/'+fn for fn in next(os.walk(fit_path))[2]]
    res_fid=open(solve_fname,'w')
    res_fid.write("# filename.fit az_c alt_c az0 alt0 a[0] a[1] a[2] b[0] b[1] b[2] c[0] c[1] c[2] d[0] d[1] d[2]\n")
    for i in range(len(fit_filenames)):
        ret = keo_solve_field_altaz(fit_filenames[i],keo_solve_pars)
        if type(ret)!=int:
            az0=ret[0]
            alt0=ret[1]
            az_c=ret[2]
            alt_c=ret[3]
            a=ret[4]
            b=ret[5]
            c=ret[6]
            d=ret[7]
            str_to_file=fit_filenames[i].split("/")[-1]+ " " + str(az_c) + " " +str(alt_c)+ " " +str(az0)+ " " +str(alt0)
            str_to_file+=" " + str(a[0]) + " " + str(a[1]) + " " + str(a[2])
            str_to_file+=" " + str(b[0]) + " " + str(b[1]) + " " + str(b[2])
            str_to_file+=" " + str(c[0]) + " " + str(c[1]) + " " + str(c[2])
            str_to_file+=" " + str(d[0]) + " " + str(d[1]) + " " + str(d[2]) + "\n"
            res_fid.write(str_to_file)
        else:
            print(ret)
            print(str(i)," ",fit_filenames[i].split("/")[-1]," ERROR")
    res_fid.close()

fit_path="../data/140824/keo"
solve_fname="solve_field_altaz_140824_keo.dat"
save_fname="keo_140824_solve.pars"
save_solve_data(fit_path,solve_fname)
save_med_solve_pars(solve_fname,save_fname)
table_fname="keo_stars_140824.table"
save_manual_fname="keo_140824_solve_manual.pars"
save_keo_manual_solve_pars(table_fname, save_manual_fname)

fit_path="../data/140826/keo"
solve_fname="solve_field_altaz_140826_keo.dat"
save_fname="keo_140826_solve.pars"
save_solve_data(fit_path,solve_fname)
save_med_solve_pars(solve_fname,save_fname)
table_fname="keo_stars_140826.table"
save_manual_fname="keo_140826_solve_manual.pars"
save_keo_manual_solve_pars(table_fname, save_manual_fname)

fit_path="../data/160829/keo"
solve_fname="solve_field_altaz_160829_keo.dat"
save_fname="keo_160829_solve.pars"
save_solve_data(fit_path,solve_fname)
save_med_solve_pars(solve_fname,save_fname)
table_fname="keo_stars_160829.table"
save_manual_fname="keo_160829_solve_manual.pars"
save_keo_manual_solve_pars(table_fname, save_manual_fname)