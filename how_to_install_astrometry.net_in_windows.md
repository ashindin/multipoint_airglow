# How to install astrometry.net under cygwin in Windows

The following instruction based on [this](https://sites.google.com/site/jmastronomy/Software/astrometry-net-setup) guide. Several things are updated.

1. **Download Cygwin ([x86](https://www.cygwin.com/setup-x86.exe) or [x86_64](https://www.cygwin.com/setup-x86.exe)). Tested with version 2.877**

2. **Install Cygwin to default directory ("C:\cygwin"). Also the following packages have to be installed with all dependencies:**  
Archive/bzip2  
Devel/gcc-core  
Devel/clang  
Devel/make  
Devel/swig  
Devel/pkg-config  
Graphics/libnetpbm-devel  
Graphics/netpbm  
Libs/libcairo-devel  
Libs/libjpeg-devel  
Libs/libgsl-devel  
Libs/libpcre-devel  
Libs/libpixman1-devel  
Libs/libX11-xcb-devel  
Libs/libxcb-glx-devel  
Libs/libXdamage-devel  
Libs/zlib-devel  
Python/python2-numpy  
Python/python2-devel  
Science/gsl  
Web/lynx  
Web/wget  
**Command to install all these packages in windows command line:**  
setup-x86.exe -P bzip2, gcc-core, clang, make, swig, pkg-config, libnetpbm-devel, netpbm, libcairo-devel, libjpeg-devel, libgsl-devel, libpcre-devel, libpixman1-devel, libX11-xcb-devel, libxcb-glx-devel, libXdamage-devel, zlib-devel, python2-numpy, python2-devel, gsl, lynx, wget  
*Note that you need to change setup-x86.exe to setup-x86_64.exe if you use the 64 bit version of cygwin*

3. **Add "C:\cygwin\bin" to PATH** 

4. **Install pip under cygwin. In cygwin terminal (as Administrator!) run:**  
lynx -source https://bootstrap.pypa.io/get-pip.py | python

5. **Install astropy via cygwin's pip:**  
CC=gcc  
pip install --no-deps astropy

6. **Install CFITSIO:**  
wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio_latest.tar.gz  
tar -xvzf cfitsio_latest.tar.gz  
cd cfitsio
./configure --prefix=/usr/local  
make  
make install  
cd ..

7. **Install astrometry.net:**  
wget http://astrometry.net/downloads/astrometry.net-latest.tar.gz  
tar -xvzf astrometry.net-latest.tar.gz  
cd astrometry.net-0.70  
cd util  
echo NETPBM_INC ?= -I /usr/include/netpbm > makefile.netpbm  
echo NETPBM_LIB ?= -L. -lnetpbm >> makefile.netpbm  
cd ..  
make CFITS_INC="-I/usr/local/include" CFITS_LIB="-L/usr/local/lib -lcfitsio"  
make extra  
make install

8. **Add astrometry binaries to cygwin's PATH:**  
echo PATH='$PATH':/usr/local/astrometry/bin >> .bash_profile