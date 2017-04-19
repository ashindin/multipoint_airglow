import numpy as np
from astropy.io import fits
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

filelists=["keo_140824_dark.files", "keo_140826_dark.files", 
"s1c_140824_dark.files", "s1c_140826_dark.files", 
"keo_160829_dark.files", "sbig_160829_dark.files"]

master_dark_fnames=[filelist[0:-10]+"masterdark.fit" for filelist in filelists]

for i in range(len(filelists)):
    filelist=filelists[i]
    fid=open(filelist,'r')
    lines=fid.readlines()[1::]
    fid.close()
    for j in range(len(lines)):
        if lines[j][-1]=='\n':
            lines[j]=lines[j][0:-1]


    hdulist = fits.open(lines[0],ignore_missing_end=True)
    img=hdulist[0].data
    hdulist.close()
    img_size=img.shape

    darks=np.zeros((img_size[0],img_size[1],len(lines)),dtype=np.int16)

    for j in range(len(lines)):
        hdulist = fits.open(lines[j],ignore_missing_end=True)
        darks[:,:,j]=hdulist[0].data
        hdulist.close()

    master_dark=np.median(darks,axis=2).astype(np.int16)

    hdu = fits.PrimaryHDU(master_dark)
    hdu.writeto(master_dark_fnames[i],overwrite=True)

