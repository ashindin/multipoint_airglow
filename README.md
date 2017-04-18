# multipoint_airglow

## Dependencies  

1. Python 3.6 + numpy, scipy, astropy, sympy and matplotlib (it may be Anaconda 4.3 or greater).
2. Ffmpeg binaries in PATH.  
3. astrometry.net 0.70 or greater (for Windows - Cygwin + astrometry.net)   

## Steps to reproduce results  

Here and below it is assumed that command "python" means python version 3.6. If your default python is 2.7 you have to change all "python" commands below to "python3". Also it is assumed that your current directory is root directory of project repository (it contains README.md file).  

### Downloading phase

1. Download and unpack source data archives: `python download.py`.  
2. In Windows and if you used git to clone this project: you have to run `python crlf.py` to correct endings in some files.  

### Astrometric calibration (solving) phase  

1. Change directory: `cd astrometric_calibration`.  
2. Solve source images (in horizontal coordinats): `python solve_field_altaz_s1c.py`, `python solve_field_altaz_keo.py` and `python solve_field_altaz_sbig.py`.  
3. Run `python solve_cam_direction_plot.py` to plot central pixel directions of solved images.
4. Examine solving results (make movie files): `python exam_solve_pars_s1c.py`, `python exam_solve_pars_keo.py`, `python exam_solve_pars_sbig.py`.  

### Spectrophotometric calibration phase  
