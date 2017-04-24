import numpy as np

def get_solve_pars(save_fname):
    a=np.zeros(3)
    b=np.zeros(3)
    c=np.zeros(3)
    d=np.zeros(3)
    fid=open(save_fname,'r')
    lines=fid.readlines()
    pars_list=lines[1].split(' ')
    pars=[float(par_str) for par_str in pars_list]
    az0,alt0,a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2],d[0],d[1],d[2]=pars
    fid.close()
    return az0,alt0,a,b,c,d

def get_scale_and_orientation_info(c,d):
    Mx=np.sqrt(c[1]**2+d[1]**2)
    My=np.sqrt(c[2]**2+d[2]**2)
#     print("Mx = ",Mx," pix/deg; My = ",My, " pix/deg")
#     print("M = ",(Mx+My)/2, " pix/deg")
#     theta=np.arctan(c[2]/d[2])
#     print("Orientation: theta = ",theta*180/np.pi, " deg")
#     gamma=np.arctan((c[1]*c[2]+d[1]*d[2])/(c[1]*d[2]-d[1]*c[2]))
#     print("Oblique: gamma = ", gamma*180/np.pi, " deg")
    return (Mx+My)/2 #, theta, gamma

def arc_get_pixel_solid_angle(M,theta):
    # theta is a Zenith angle (theta==0 at Zenith)
    A=M/np.pi*180
    if theta==0. or theta==0:
        Omega=1/A**2
    else:
        Omega=1/A**2*(np.sin(theta)/theta)
    return Omega
def tan_get_pixel_solid_angle(M,theta):
    # theta is a Zenith angle (theta==0 at Zenith)
    A=M/np.pi*180
    Omega=1/A**2*(np.cos(theta)**3)
    return Omega

def get_brightness_in_Rayleighs(dlambda,sol_angle,num_of_pixels,flux):
    const=4*np.pi*6.3*dlambda/6.626/3/(sol_angle*10**6)
#     print(const)
    R=(const*10**4)*flux/num_of_pixels
    return R
    
# az0,alt0,a,b,c,d = get_solve_pars("../astrometric_calibration/keo_140824_solve_manual.pars")
# M_keo = get_scale_and_orientation_info(c,d)
# sol_angle_keo=arc_get_pixel_solid_angle(M_keo,36*np.pi/180)
# print(get_brightness_in_Rayleighs(M_keo,20,sol_angle_keo,6,236.5*10**-3))

# az0,alt0,a,b,c,d = get_solve_pars("../astrometric_calibration/s1c_140824_solve.pars")
# M_s1c = get_scale_and_orientation_info(c,d)
# sol_angle_s1c=tan_get_pixel_solid_angle(M_s1c,5*np.pi/180)
# get_brightness_in_Rayleighs(M_s1c,100,sol_angle_s1c,12,236.5*10**-3)

def load_sp_catalog(cat_fname):
    fid=open(cat_fname,'r')
    lines=fid.readlines()[2::]
    fid.close()
    num_lines=len(lines)
    NUM=np.zeros(num_lines,np.int)
    BS_ID=np.zeros(num_lines,np.int)
    RA=np.zeros(num_lines)
    DEC=np.zeros(num_lines)
    MAG=np.zeros(num_lines)
    FLUX=np.zeros(num_lines)
    SP_type=[]
    for i in range(num_lines):
        line_list=lines[i].split(" ")
        NUM[i]=int(line_list[0])
        BS_ID[i]=int(line_list[1])
        RA[i]=float(line_list[2])*np.pi/180
        DEC[i]=float(line_list[3])*np.pi/180
        MAG[i]=float(line_list[4])
        FLUX[i]=float(line_list[5])
        SP_type.append(line_list[6][0:-1])
    return NUM, BS_ID, RA, DEC, MAG, FLUX, SP_type