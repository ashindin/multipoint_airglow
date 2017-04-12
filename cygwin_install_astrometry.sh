#!/bin/bash
cd ~
lynx -source https://bootstrap.pypa.io/get-pip.py | python

CC=gcc  
pip install --no-deps astropy

wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio_latest.tar.gz  
tar -xvzf cfitsio_latest.tar.gz  
cd cfitsio
./configure --prefix=/usr/local  
make  
make install  
cd .

tar -xvzf astrometry.net-latest.tar.gz  
cd astrometry.net-0.70  
cd util  
echo NETPBM_INC ?= -I /usr/include/netpbm > makefile.netpbm  
echo NETPBM_LIB ?= -L. -lnetpbm >> makefile.netpbm  
cd ..  
make CFITS_INC="-I/usr/local/include" CFITS_LIB="-L/usr/local/lib -lcfitsio"  
make extra  
make install

echo PATH='$PATH':/usr/local/astrometry/bin >> .bash_profile