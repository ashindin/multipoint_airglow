
# coding: utf-8

# In[1]:

import datetime
import os, sys, inspect
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.signal as ss

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"../astrometric_calibration")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from tan_module import *
from sp_cal_module import *


# In[2]:

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
def make_movie_from_pngs(png_prefix, num_frames, movie_fname):
    com_line="ffmpeg -y -r 5 -f image2 -s 1280x720 -i " + png_prefix + "%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p " +"-vframes "+ str(num_frames) +" "+ movie_fname
#     print(com_line)
    os.system(com_line)
    return None


# In[3]:

cat_fname="SP_CATALOG.csv"
NUM, BS_ID, RA, DEC, MAG, FLUX, SP_type = load_sp_catalog(cat_fname)
dx=0.05
ms=10
mew=1


# In[4]:

np.set_printoptions(precision=1,linewidth=100,suppress=True)


# In[11]:

def s1c_sp_calibration(fit_path,masterdark_fname,masterflat_fname,solve_pars_fname,save_fname=None,area_rad=4, med_size=31,lat_deg=56.1501667,lon_deg=46.1050833,hei_m=183.):
    if save_fname==None:
        save_fname=solve_pars_fname.split('/')[-1][0:-11]+'.spcal'

    hdulist = fits.open(masterdark_fname,ignore_missing_end=True)
    masterdark=hdulist[0].data
    hdulist.close()

    hdulist = fits.open(masterflat_fname,ignore_missing_end=True)
    masterflat=hdulist[0].data
    hdulist.close()

    err,az0,alt0,a,b,c,d=get_solve_pars(solve_pars_fname)
    M_s1c = get_scale_and_orientation_info(c,d)

    spath="./.temp/"
    if not os.path.exists(spath):
        os.makedirs(spath)
    fit_filenames=[fit_path+'/'+fn for fn in next(os.walk(fit_path))[2]]

    R_median=np.zeros(len(fit_filenames))
    R_std=np.zeros(len(fit_filenames))

    fig=plt.figure(figsize=(12.8,7.2))
    fig.set_size_inches(12.8, 7.2)

#     ax=plt.axes(position=[0.035000000000000003+dx/2, 0.061764705882352888, 0.50624999999999998, 0.89999999999999991])
    ax=plt.axes(position=[0.035000000000000003+dx/2, 0.26, 0.50624999999999998, 0.7])
    plt.axis('off')
    ax1=plt.axes(position=[0.58205882352941185+dx, 0.71052941176470596-0.15, 0.35,0.4])
    plt.grid()
    ax2=plt.axes(position=[0.58205882352941185+dx, 0.061764705882352888, 0.35, 0.4])
    plt.grid()

    ax3=plt.axes(position=[0.035000000000000003+dx/2,0.04, 0.5/4, 0.7/4])
    plt.grid(b=True)
    ax4=plt.axes(position=[0.035000000000000003+dx/2+0.19,0.04, 0.5/4, 0.7/4])
    plt.grid(b=True)
    ax5=plt.axes(position=[0.035000000000000003+dx/2+0.38,0.04, 0.5/4, 0.7/4])
    plt.grid(b=True)

    R_day=0

    vrange=1000
    vshift=2000

    alt_min=83*np.pi/180
#     area_rad=3
#     med_size=31

    XX,YY=np.meshgrid(range(2*area_rad+1),range(2*area_rad+1))

    fid=open(save_fname,'w')
    fid.write("# Camera calibration coefficients [Rayleighs per ADC unit] for each fit file:\n")

    for i in range(len(fit_filenames)):
#     for i in range(5):
#     for i in range(58,59):
        sys.stdout.write('\r')
        sys.stdout.write("Processing frame "+str(i+1)+"/"+str(len(fit_filenames)))
        sys.stdout.flush()
        fit_fname=fit_filenames[i]
        hdulist = fits.open(fit_fname,ignore_missing_end=True)
        img=hdulist[0].data.astype('float')
        img0=np.copy(img)
        img=(img-masterdark.astype('float'))/masterflat.astype('float')
        img_medfilt=img - ss.medfilt(img,kernel_size=med_size)
#         np.save('img.npy',img)
        hdulist.close()

        med_img=np.median(img)
        max_img0=np.max(img0)
#         print(max_img0)

        s1c_site=EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg, height=hei_m*u.m)
        date_obs=s1c_get_date_obs(fit_fname,ut_shift=-4)

        png_prefix=spath+"frame_s1c_"+str(date_obs.year)[2::]+"{0:0>2}".format(date_obs.month)+"{0:0>2}".format(date_obs.day)+"_"

        BS_coord=SkyCoord(RA, DEC, frame='icrs', unit='rad');
        altaz=BS_coord.transform_to(AltAz(obstime=date_obs, location=s1c_site,temperature=20*u.deg_C,pressure=1013*u.hPa,
                                               relative_humidity=0.5,obswl=630.0*u.nm))
        AZ=altaz.az.rad
        ALT=altaz.alt.rad

        X,Y = tan_hor2pix(AZ,ALT,az0,alt0,c,d)

        # Catalog filtration
        filt_mask=np.zeros(len(NUM),dtype=bool)
        for j in range(len(NUM)):
            if ALT[j]>=alt_min and X[j]>=10 and Y[j]>=10 and X[j]<=img.shape[1]-10 and Y[j]<=img.shape[0]-10:
                filt_mask[j]=True
        NUM_filt=NUM[filt_mask]
        BS_ID_filt=BS_ID[filt_mask]
        RA_filt=RA[filt_mask]
        DEC_filt=DEC[filt_mask]
        MAG_filt=MAG[filt_mask]
        FLUX_filt=FLUX[filt_mask]
        SP_type_filt=[SP_type[j] for j in range(len(SP_type)) if filt_mask[j]==True]
        ALT_filt=ALT[filt_mask]
        AZ_filt=AZ[filt_mask]
        X_filt=X[filt_mask]
        Y_filt=Y[filt_mask]

        star_pixels=np.zeros(len(NUM_filt),dtype=int)
        star_adc=np.zeros(len(NUM_filt))

        AREA=np.zeros((2*area_rad+1,2*area_rad+1,len(NUM_filt)))
        AREA0=np.zeros((2*area_rad+1,2*area_rad+1,len(NUM_filt)))

        for j in range(len(NUM_filt)):
            sum_temp=0
            num=0
            AREA[:,:,j]=np.copy(img_medfilt[int(Y_filt[j])-area_rad:int(Y_filt[j])+area_rad+1, int(X_filt[j])-area_rad:int(X_filt[j])+area_rad+1])
            AREA0[:,:,j]=np.copy(img0[int(Y_filt[j])-area_rad:int(Y_filt[j])+area_rad+1, int(X_filt[j])-area_rad:int(X_filt[j])+area_rad+1])

            arg_max=np.argmax(AREA[:,:,j])
#             print("max=", np.max(AREA[:,:,j]))
            x0=XX.flat[arg_max]
            y0=YY.flat[arg_max]

            AREA[:,:,j]=np.copy(img_medfilt[int(Y_filt[j])-area_rad+y0-area_rad:int(Y_filt[j])+y0+1, int(X_filt[j])-area_rad+x0-area_rad:int(X_filt[j])+x0+1])
            AREA0[:,:,j]=np.copy(img0[int(Y_filt[j])-area_rad+y0-area_rad:int(Y_filt[j])+y0+1, int(X_filt[j])-area_rad+x0-area_rad:int(X_filt[j])+x0+1])

            Rast=np.sqrt((XX-area_rad)**2+(YY-area_rad)**2)

#             print(Rast)

            area=np.copy(AREA[:,:,j])

            for k in range((2*area_rad+1)**2):
                if Rast.flat[k]<=area_rad:
                    num+=1
                    sum_temp+=AREA[:,:,j].flat[k]
            star_pixels[j]=num
            star_adc[j]=sum_temp
#             print(star_pixels[j],star_adc[j])
#             print(AREA[:,:,j])
#             print(" ")


        # Catalog filtration 2
        filt_mask=np.zeros(len(NUM_filt),dtype=bool)
        for j in range(len(NUM_filt)):
            if np.max(AREA0[:,:,j])<0.9*max_img0:
                filt_mask[j]=True

        AREA_filt2=np.copy(AREA[:,:,filt_mask])
        AREA0_filt2=np.copy(AREA0[:,:,filt_mask])

        NUM_filt2=NUM_filt[filt_mask]
        BS_ID_filt2=BS_ID_filt[filt_mask]
        RA_filt2=RA_filt[filt_mask]
        DEC_filt2=DEC_filt[filt_mask]
        MAG_filt2=MAG_filt[filt_mask]
        FLUX_filt2=FLUX_filt[filt_mask]
        SP_type_filt2=[SP_type_filt[j] for j in range(len(SP_type_filt)) if filt_mask[j]==True]
        X_filt2=X_filt[filt_mask]
        Y_filt2=Y_filt[filt_mask]
        ALT_filt2=ALT_filt[filt_mask]
        AZ_filt2=AZ_filt[filt_mask]
        star_pixels_filt2=star_pixels[filt_mask]
        star_adc_filt2=star_adc[filt_mask]

        R=np.zeros(len(NUM_filt2))
        for j in range(len(NUM_filt2)):
            sol_angle=tan_get_pixel_solid_angle(M_s1c,np.pi/2-ALT_filt2[j])
            br=get_brightness_in_Rayleighs(100,sol_angle,1,FLUX_filt2[j])
            sa=star_adc_filt2[j]
            R[j]=br/sa
#             print('br[R]=', br, "; sa[ADC.u]=",sa,"; R=",R[j])

        x_bord=np.arange(2*area_rad+1,dtype=float)
        y1_bord=np.sqrt(area_rad**2-(x_bord-area_rad)**2)+area_rad+0.5
        y2_bord=-np.sqrt(area_rad**2-(x_bord-area_rad)**2)+area_rad+0.5

        if len(NUM_filt2)>0:

            sort_ord=np.argsort(star_pixels_filt2)

            plt.sca(ax3)
            idd=0
            area1=AREA_filt2[:,:,sort_ord[idd]]
            area1_max=np.max(AREA0_filt2[:,:,sort_ord[idd]])
            area1_id=BS_ID_filt2[sort_ord[idd]]
            plt.pcolormesh(area1, cmap="seismic",vmin=-1000, vmax=1000)

            plt.plot(x_bord+0.5,y1_bord,'k-',lw=2)
            plt.plot(x_bord+0.5,y2_bord,'k-',lw=2)

            plt.ylim(area_rad*2+1,0)
            plt.xlim(0,area_rad*2+1)
            plt.title(str(area1_id),loc='left',fontsize='smaller')
            plt.title(str(int(star_adc_filt2[sort_ord[idd]])),loc='right',fontsize='smaller')
            plt.title(str(int(R[sort_ord[idd]]*1000)/1000),fontsize='smaller')
            plt.ylabel(str(int(R[sort_ord[idd]]*star_adc_filt2[sort_ord[idd]])))

            if len(NUM_filt2)>1:
                plt.sca(ax4)
                idd=1
                area2=AREA_filt2[:,:,sort_ord[idd]]
                area2_max=np.max(AREA0_filt2[:,:,sort_ord[idd]])
                area2_id=BS_ID_filt2[sort_ord[idd]]
                plt.pcolormesh(area2, cmap="seismic",vmin=-1000, vmax=1000)

                plt.plot(x_bord+0.5,y1_bord,'k-',lw=2)
                plt.plot(x_bord+0.5,y2_bord,'k-',lw=2)

                plt.ylim(area_rad*2+1,0)
                plt.xlim(0,area_rad*2+1)
                plt.title(str(area2_id),loc='left',fontsize='smaller')
                plt.title(str(int(star_adc_filt2[sort_ord[idd]])),loc='right',fontsize='smaller')
                plt.title(str(int(R[sort_ord[idd]]*1000)/1000),fontsize='smaller')
                plt.ylabel(str(int(R[sort_ord[idd]]*star_adc_filt2[sort_ord[idd]])))
            if len(NUM_filt2)>2:
                plt.sca(ax5)
                idd=2
                area3=AREA_filt2[:,:,sort_ord[idd]]
                area3_max=np.max(AREA0_filt2[:,:,sort_ord[idd]])
                area3_id=BS_ID_filt2[sort_ord[idd]]
                plt.pcolormesh(area3, cmap="seismic",vmin=-1000, vmax=1000)

                plt.plot(x_bord+0.5,y1_bord,'k-',lw=2)
                plt.plot(x_bord+0.5,y2_bord,'k-',lw=2)

                plt.ylim(area_rad*2+1,0)
                plt.xlim(0,area_rad*2+1)
                plt.title(str(area3_id),loc='left',fontsize='smaller')
                plt.title(str(int(star_adc_filt2[sort_ord[idd]])),loc='right',fontsize='smaller')
                plt.title(str(int(R[sort_ord[idd]]*1000)/1000),fontsize='smaller')
                plt.ylabel(str(int(R[sort_ord[idd]]*star_adc_filt2[sort_ord[idd]])))



        # print(len(R),R)
        if len(R)>0:
            R_median[i]=np.median(R)
            temp=(R-R_median[i])**2
            R_std[i]=np.sqrt(np.median(temp))
        #     print(R_median[i])

        plt.sca(ax1)
        plt.plot([0, 50000], [R_median[i], R_median[i]], c='r', lw=2)
        plt.scatter(R*star_adc_filt2,R,s=np.pi*(star_adc_filt2/12000*7)**2)
#         print(R[0]*star_adc_filt2[0])
#         print(R[0]*star_adc_filt2[0]*star_pixels_filt2[0])

        for j in range(len(X_filt2)):
            plt.text(R[j]*star_adc_filt2[j],R[j],str(BS_ID_filt2[j])+"_"+str(star_pixels_filt2[j]))
        plt.ylabel('Calibration coef. [R/ADCu]')
        plt.xlabel('Star summ brightness [R]')
        plt.title(date_obs, loc='left')
        plt.title(R_median[i], loc='right')
        ax1.set_xlim((0, 30000))
        ax1.set_ylim((0, 5))
        plt.grid(b=True)

        plt.sca(ax2)
        plt.xlim([0,len(fit_filenames)-1])
        plt.ylim([0,5])
        plt.plot(i,R_median[i],"r.")
        plt.grid(b=True)
        plt.ylabel('Calibration coef. [R/ADCu]')
        plt.xlabel('frame number')
        if i==len(fit_filenames)-1:
            R_day=np.median(R_median[np.where(R_median>0)])
            plt.plot([0, 5000], [R_day, R_day], c='b', lw=2)
            plt.title(R_day, loc='right')

        plt.sca(ax)
        plt.pcolormesh(img, cmap="gray", vmin=np.median(img)-100, vmax=np.median(img)+100)
        plt.axis('equal')
        plt.plot(X,Y,marker="o", lw=0.,mew=mew,mec="b", mfc='none',ms=ms)
        plt.plot(X_filt,Y_filt,marker="o", lw=0.,mew=mew,mec="g", mfc='none',ms=ms)
        plt.plot(X_filt2,Y_filt2,marker="o", lw=0.,mew=mew,mec="r", mfc='none',ms=ms)
        for j in range(len(X_filt2)):
            plt.text(X_filt2[j],Y_filt2[j],str(BS_ID_filt2[j])+"_"+ "{0:4.2f}".format(R[j]) ,color='w')
        ax.set_xlim((1,372))
        ax.set_ylim((281,1))
        plt.title(fit_fname.split('/')[-1] + " " + "{0:0>2}".format(date_obs.hour) + ":" + "{0:0>2}".format(date_obs.minute) + ":" + "{0:0>2}".format(date_obs.second))
    #     plt.show()
        plt.axis('off')
        png_fname=png_prefix+"{0:0>4}".format(i+1)+".png"
    #   print(png_fname)
        plt.savefig(png_fname)
        ax.clear()
        ax1.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()

        fid.write(fit_fname.split('/')[-1] + " " + str(R_median[i])+ " " + str(R_std[i]) + "\n")

    plt.close()
    sys.stdout.write('\n')
    sys.stdout.flush()
    fid.close()

    return png_prefix, len(fit_filenames), R_median, R_std


# In[12]:

fit_path="../data/140824/s1c"
movie_fname="s1c_140824_calibration.mp4"
solve_pars_fname="../astrometric_calibration/s1c_140824_solve.pars"
masterdark_fname="s1c_140824_masterdark.fit"
masterflat_fname="s1c_master.flat"
save_fname="s1c_140824_day.spcal"
png_prefix, num_frames, R_median, R_std = s1c_sp_calibration(fit_path,masterdark_fname,masterflat_fname,solve_pars_fname,area_rad=4,med_size=21)
make_movie_from_pngs(png_prefix, num_frames, movie_fname)
R_filt=R_median[np.where(R_median>0)]
R_day=np.median(R_filt)
R_std2=np.std(R_filt)
print(R_day, R_std2)
fid=open(save_fname,'w')
fid.write("# Median camera calibration coefficient [Rayleighs per ADC unit] its std for 14/08/24:\n")
fid.write(str(R_day)+" "+str(R_std2))
fid.close()


# In[13]:

fit_path="../data/140826/s1c"
movie_fname="s1c_140826_calibration.mp4"
solve_pars_fname="../astrometric_calibration/s1c_140826_solve.pars"
masterdark_fname="s1c_140826_masterdark.fit"
masterflat_fname="s1c_master.flat"
save_fname="s1c_140826_day.spcal"
png_prefix, num_frames, R_median, R_std = s1c_sp_calibration(fit_path,masterdark_fname,masterflat_fname,solve_pars_fname,area_rad=4,med_size=21)
make_movie_from_pngs(png_prefix, num_frames, movie_fname)
R_median=R_median[500::]
R_filt=R_median[np.where(R_median>0)]
R_day=np.median(R_filt)
R_std2=np.std(R_filt)
print(R_day, R_std2)
fid=open(save_fname,'w')
fid.write("# Median camera calibration coefficient [Rayleighs per ADC unit] and its std for 14/08/26:\n")
fid.write(str(R_day)+" "+str(R_std2))
fid.close()

