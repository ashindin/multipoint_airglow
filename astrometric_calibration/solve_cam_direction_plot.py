
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from save_solve_pars import *


# In[3]:

solve_fnames=[
    "solve_field_altaz_140824_s1c.dat",
    "solve_field_altaz_140826_s1c.dat",
    "solve_field_altaz_160829_sbig.dat",
    "solve_field_altaz_140824_keo.dat",
    "solve_field_altaz_140826_keo.dat",
    "solve_field_altaz_160829_keo.dat",
]


# In[10]:

for solve_fname in solve_fnames:
    az_c,alt_c,az0,alt0,a,b,c,d=get_solve_coefs(solve_fname)
    med_ind= get_med_ind(az_c,alt_c)
    fig = plt.figure(figsize=(9,9))
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_zero_location("N")
    plt.plot(az_c,(np.pi/2-alt_c)*180/np.pi,'bo')
    plt.plot(az_c[med_ind],(np.pi/2-alt_c[med_ind])*180/np.pi,'ro')
    plt.title(solve_fname)
    if solve_fname=="solve_field_altaz_160829_sbig.dat": ax.set_rmax(3)
#     plt.show()
    plt.savefig(solve_fname[:-3]+'png')
    plt.close()


# In[ ]:



