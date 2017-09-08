import datetime
import numpy as np
import os
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import scipy.optimize as so

def get_solve_coefs(solve_fname):
    fid=open(solve_fname,'r')
    lines=fid.readlines()
    fid.close()
    lines=lines[1::]
    num_lines=len(lines)
    fit_names=[]
    az_c=np.zeros(num_lines)
    alt_c=np.zeros(num_lines)
    az0=np.zeros(num_lines)
    alt0=np.zeros(num_lines)
    err=np.zeros(num_lines)
    a=np.zeros((num_lines,3))
    b=np.zeros((num_lines,3))
    c=np.zeros((num_lines,3))
    d=np.zeros((num_lines,3))
    for i in range(0,num_lines):
        line_list=lines[i].split(' ')
        fit_names.append(line_list[0])
        err[i]=float(line_list[1])
        az_c[i]=float(line_list[2])
        alt_c[i]=float(line_list[3])        
        az0[i]=float(line_list[4])
        alt0[i]=float(line_list[5])        
        a[i,0]=float(line_list[6])
        a[i,1]=float(line_list[7])
        a[i,2]=float(line_list[8])
        b[i,0]=float(line_list[9])
        b[i,1]=float(line_list[10])
        b[i,2]=float(line_list[11])
        c[i,0]=float(line_list[12])
        c[i,1]=float(line_list[13])
        c[i,2]=float(line_list[14])
        d[i,0]=float(line_list[15])
        d[i,1]=float(line_list[16])
        d[i,2]=float(line_list[17])
    return err,az_c,alt_c,az0,alt0,a,b,c,d
    
def get_med_ind(az_c,alt_c):
    x= np.cos(alt_c)*np.cos(az_c)
    y= np.cos(alt_c)*np.sin(az_c)
    x_med=np.median(x)
    y_med=np.median(y)
    r=np.sqrt((x-x_med)**2+(y-y_med)**2)
#     print(r)
    return np.argmin(r)

def save_med_solve_pars(solve_fname,save_fname):    
    err,az_c,alt_c,az0,alt0,a,b,c,d=get_solve_coefs(solve_fname)
    med_ind= get_med_ind(az_c,alt_c)
    err_med=err[med_ind]
    az0_med=az0[med_ind]
    alt0_med=alt0[med_ind]
    a_med=a[med_ind,:]
    b_med=b[med_ind,:]
    c_med=c[med_ind,:]
    d_med=d[med_ind,:]
    fid=open(save_fname,'w')
    fid.write("# error_pix az0 alt0 a[0] a[1] a[2] b[0] b[1] b[2] c[0] c[1] c[2] d[0] d[1] d[2]\n")
    str_to_file=str(err_med)+" "+str(az0_med)+ " " +str(alt0_med)
    str_to_file+=" " + str(a_med[0]) + " " + str(a_med[1]) + " " + str(a_med[2])
    str_to_file+=" " + str(b_med[0]) + " " + str(b_med[1]) + " " + str(b_med[2])
    str_to_file+=" " + str(c_med[0]) + " " + str(c_med[1]) + " " + str(c_med[2])
    str_to_file+=" " + str(d_med[0]) + " " + str(d_med[1]) + " " + str(d_med[2]) + "\n"
    fid.write(str_to_file)
    fid.close()

def keo_get_date_obs(filename):
    fid_fit=fits.open(filename);
    date_obs_str=fid_fit[0].header["DATE-OBS"]  
    fid_fit.close()
    return datetime.datetime.strptime(date_obs_str,"%Y-%m-%dT%H:%M:%S")

import arc_module

def save_keo_manual_solve_pars(stars_fname,save_fname,lat_deg=55.9305361,lon_deg=48.7444861,hei_m=91.):
    
    keo_site=EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg, height=hei_m*u.m)
    
    fid=open(stars_fname,'r')
    fit_fname=fid.readline()[2:-1]
    fid.close() 
    
    date_obs=keo_get_date_obs(fit_fname)
    
    table=np.loadtxt(stars_fname)
    
    X=table[:,0]
    Y=table[:,1]
    sao_name=[]
    
    AZ=np.zeros(table.shape[0])
    ALT=np.zeros(table.shape[0])
    
    for i in range(table.shape[0]):
        sao_name="SAO " + str(int(table[i,2]))
        sc=SkyCoord.from_name(sao_name)
        aa=sc.transform_to(AltAz(obstime=date_obs, location=keo_site,temperature=20*u.deg_C,pressure=1013*u.hPa,
                                               relative_humidity=0.5,obswl=630.0*u.nm))
        AZ[i]=aa.az.rad
        ALT[i]=aa.alt.rad
    
#     print(np.vstack((table[:,2].T,X,Y,AZ*180/np.pi,ALT*180/np.pi)).T)
    
    res=so.minimize(arc_module.arc_calc_pix2st_coefs_discrep,(0,np.pi/2),(AZ,ALT,X,Y),method='Nelder-Mead')
    res_fun=res.fun
    res_x=res.x
    if res_x[1]>np.pi/2:
        res_x[1]=np.pi-res_x[1]
        res_x[0]=res_x[0]+np.pi
    if res_x[0]<0:
        res_x[0]=res_x[0]+2*np.pi
    if res_x[0]>2*np.pi:
        res_x[0]=res_x[0]%(2*np.pi)

    az0_med=res_x[0];
    alt0_med=res_x[1];
    
    a_med,b_med,c_med,d_med=arc_module.arc_calc_pix2st_coefs((az0_med,alt0_med),AZ,ALT,X,Y)
    #print("sigma = ", arc_module.arc_calc_pix2st_coefs_discrep((az0,alt0),AZ,ALT,X,Y))
    
    fid=open(save_fname,'w')
    fid.write("# error_pix az0 alt0 a[0] a[1] a[2] b[0] b[1] b[2] c[0] c[1] c[2] d[0] d[1] d[2]\n")
    str_to_file=str(res_fun) + " " + str(az0_med)+ " " +str(alt0_med)
    str_to_file+=" " + str(a_med[0]) + " " + str(a_med[1]) + " " + str(a_med[2])
    str_to_file+=" " + str(b_med[0]) + " " + str(b_med[1]) + " " + str(b_med[2])
    str_to_file+=" " + str(c_med[0]) + " " + str(c_med[1]) + " " + str(c_med[2])
    str_to_file+=" " + str(d_med[0]) + " " + str(d_med[1]) + " " + str(d_med[2]) + "\n"
    fid.write(str_to_file)
    fid.close()
    
    #x_new,y_new = arc_module.arc_hor2pix(AZ,ALT,az0,alt0,c,d)
    #print(np.vstack((table[:,2].T,X,Y,AZ*180/np.pi,ALT*180/np.pi,x_new,y_new)).T)
    
    az_c, alt_c = arc_module.arc_pix2hor(255,255,az0_med,alt0_med,a_med,b_med)
    print("Central pixel direction (AZ, ALT in deg):",az_c*180/np.pi,alt_c*180/np.pi)       
    #return az0,alt0,az_c,alt_c,a,b,c,d
