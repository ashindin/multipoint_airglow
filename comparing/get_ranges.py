import numpy as np
import scipy.optimize as so

def read_ellipsoid_args(dat_fn):
    fid=open(dat_fn,'r')
    fid.readline()
    success_str=fid.readline()
    status_str=fid.readline()
    #~ if success_str!="SUCCESS = True\n" or status_str!="STATUS = 0\n":
        #~ fid.close()
        #~ return None
    #~ else:
    fid.readline()
    fid.readline()
    fid.readline()
    args_list=fid.readline().split(' ')
    return float(args_list[4]), float(args_list[5]), float(args_list[6])


def read_spheroid_incl_args(dat_fn):
    fid=open(dat_fn,'r')
    fid.readline()
    success_str=fid.readline()
    status_str=fid.readline()
    #~ if success_str!="SUCCESS = True\n" or status_str!="STATUS = 0\n":
        #~ fid.close()
        #~ return None
    #~ else:
    fid.readline()
    fid.readline()
    fid.readline()
    args_list=fid.readline().split(' ')
    return float(args_list[4]), float(args_list[5]), float(args_list[6]), float(args_list[7])
    
def get_spheroid_incl_range(dat_fn):
    
    def z_surf1(arg, a2, a3, a4, a5):
        x=arg[0]
        y=arg[1]    
        surf=(-a3*np.sqrt(-a2**4*np.sin(a5)**2 + a2**4 + a2**2*a3**2*np.sin(a5)**2 + a2**2*x**2*np.sin(a4)**2*np.sin(a5)**2 - a2**2*x**2 - 2*a2**2*x*y*np.sin(a4)*np.sin(a5)**2*np.cos(a4) - a2**2*y**2*np.sin(a4)**2*np.sin(a5)**2 + a2**2*y**2*np.sin(a5)**2 - a2**2*y**2 - a3**2*x**2*np.sin(a4)**2*np.sin(a5)**2 + 2*a3**2*x*y*np.sin(a4)*np.sin(a5)**2*np.cos(a4) + a3**2*y**2*np.sin(a4)**2*np.sin(a5)**2 - a3**2*y**2*np.sin(a5)**2) + (a2 - a3)*(a2 + a3)*(x*np.cos(a4) + y*np.sin(a4))*np.sin(a5)*np.cos(a5))/(a2**2*np.cos(a5)**2 + a3**2*np.sin(a5)**2)
        return surf
        
    def z_surf2(arg, a2, a3, a4, a5):
        x=arg[0]
        y=arg[1]    
        surf=(a3*np.sqrt(-a2**4*np.sin(a5)**2 + a2**4 + a2**2*a3**2*np.sin(a5)**2 + a2**2*x**2*np.sin(a4)**2*np.sin(a5)**2 - a2**2*x**2 - 2*a2**2*x*y*np.sin(a4)*np.sin(a5)**2*np.cos(a4) - a2**2*y**2*np.sin(a4)**2*np.sin(a5)**2 + a2**2*y**2*np.sin(a5)**2 - a2**2*y**2 - a3**2*x**2*np.sin(a4)**2*np.sin(a5)**2 + 2*a3**2*x*y*np.sin(a4)*np.sin(a5)**2*np.cos(a4) + a3**2*y**2*np.sin(a4)**2*np.sin(a5)**2 - a3**2*y**2*np.sin(a5)**2) + (a2 - a3)*(a2 + a3)*(x*np.cos(a4) + y*np.sin(a4))*np.sin(a5)*np.cos(a5))/(a2**2*np.cos(a5)**2 + a3**2*np.sin(a5)**2)
        return -surf
    
    a2, a3, a4, a5 = read_spheroid_incl_args(dat_fn)
    res1 = so.minimize(z_surf1, [0.,0.], (a2, a3, a4, a5), method='Nelder-Mead')    
    dzb=abs(res1.fun)
    res2 = so.minimize(z_surf2, [0.,0.], (a2, a3, a4, a5), method='Nelder-Mead')
    dzt=abs(res2.fun)
    return dzb, dzt

def get_ellipsoid_range(dat_fn):    
    a3, a4, dzb = read_ellipsoid_args(dat_fn)
    dzt=dzb
    return dzb, dzt

#~ dat_fn='../solving/spheroid_incl/140824_175227_000000_spheroid_incl.dat'
#~ print(get_spheroid_incl_range(dat_fn))

def read_drop_args(dat_fn):
    fid=open(dat_fn,'r')
    fid.readline()
    success_str=fid.readline()
    status_str=fid.readline()
    #~ if success_str!="SUCCESS = True\n" or status_str!="STATUS = 0\n":
        #~ fid.close()
        #~ return None
    #~ else:
    fid.readline()
    fid.readline()
    fid.readline()
    args_list=fid.readline().split(' ')
    
    aa=float(args_list[7])
    
    if aa<0:
        a5=2.
        a6=6.
    if aa>1.:
        a5=6.
        a6=2.
    if aa<0.5:
        a5=2.+aa*8
        a6=6.
    else:
        a5=6.
        a6=6.-(aa*8-4.)
        
        
    return float(args_list[6]), a5, a6, aa


def get_drop_range(dat_fn):
    
    def fun_to_minimize(z, a4, a5, a6):    
        fun=(z/a4)**(a5 - 1)*((a4 - z)/a4)**(a6 - 1)*((a5 - 1)/(a5 + a6 - 2))**(-a5)*((a6 - 1)/(a5 + a6 - 2))**(-a6)*(a5*a6 - a5 - a6 + 1)/(a5**2 + 2*a5*a6 - 4*a5 + a6**2 - 4*a6 + 4) - np.exp(-1)
        return -fun
        
    def fun_to_solve(z, a4, a5, a6):    
        fun=(z/a4)**(a5 - 1)*((a4 - z)/a4)**(a6 - 1)*((a5 - 1)/(a5 + a6 - 2))**(-a5)*((a6 - 1)/(a5 + a6 - 2))**(-a6)*(a5*a6 - a5 - a6 + 1)/(a5**2 + 2*a5*a6 - 4*a5 + a6**2 - 4*a6 + 4) - np.exp(-1)
        return fun
    
    a4, a5, a6, aa = read_drop_args(dat_fn)    
    
    res=so.minimize(fun_to_minimize, a4/2, (a4, a5, a6), method='Nelder-Mead')
    max_z=res.x[0]
    
    #~ res=so.root(fun_to_solve, a4/2, (a4, a5, a6))
    res=so.root(fun_to_solve, 0.5*max_z, (a4, a5, a6))
    left_z=res.x[0]
    
    #~ res=so.root(fun_to_solve, 4/5*a4, (a4, a5, a6))
    res=so.root(fun_to_solve, 1.5*max_z, (a4, a5, a6))
    right_z=res.x[0]
    
    dzb=max_z-left_z
    dzt=right_z-max_z
    if dzb < 0 or dzt < 0:
        print(dat_fn)
        
    
    if aa>=0.5 and dzb>=dzt:    
        return dzb, dzt
    elif aa<0.5 and dzb<dzt:    
        return dzb, dzt
    else:
        print("NOT CORRECT SOLVING:", dat_fn, aa, dzb,dzt)
        return dzb, dzt

#~ dat_fn='../solving/drop/160829_174512_500500_drop.dat'
#~ print(read_drop_args(dat_fn))
#~ print(get_drop_range(dat_fn))

#~ def fun_to_minimize(z, a4, a5, a6):    
	#~ fun=(z/a4)**(a5 - 1)*((a4 - z)/a4)**(a6 - 1)*((a5 - 1)/(a5 + a6 - 2))**(-a5)*((a6 - 1)/(a5 + a6 - 2))**(-a6)*(a5*a6 - a5 - a6 + 1)/(a5**2 + 2*a5*a6 - 4*a5 + a6**2 - 4*a6 + 4) - np.exp(-1)
	#~ return -fun
#~ def fun_to_solve(z, a4, a5, a6):    
	#~ fun=(z/a4)**(a5 - 1)*((a4 - z)/a4)**(a6 - 1)*((a5 - 1)/(a5 + a6 - 2))**(-a5)*((a6 - 1)/(a5 + a6 - 2))**(-a6)*(a5*a6 - a5 - a6 + 1)/(a5**2 + 2*a5*a6 - 4*a5 + a6**2 - 4*a6 + 4) - np.exp(-1)
	#~ return fun

#~ a4, a5, a6, aa = read_drop_args(dat_fn) 
#~ z_axe=np.linspace(-100000,100000,2000)
#~ fun=fun_to_solve(z_axe,a4,a5,a6)
#~ print(z_axe,fun)

#~ res=so.minimize(fun_to_minimize, a4/2, (a4, a5, a6), method='Nelder-Mead')
#~ max_z=res.x[0]
#~ print(max_z)

#~ res=so.root(fun_to_solve, a4/2, (a4, a5, a6))
#~ left_z=res.x[0]
#~ print(left_z)

#~ import matplotlib.pyplot as plt

#~ plt.plot(z_axe-max_z,fun)
#~ plt.grid()
#~ plt.show()
