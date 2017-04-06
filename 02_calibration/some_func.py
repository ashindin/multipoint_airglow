import numpy as np

def cart2spheric(x,y,z):
    r=np.sqrt(x*x+y*y+z*z)
    phi=np.arctan2(y,x)
    theta=np.arccos(z/r)
    return (r, phi, theta)

def spheric2cart(r,phi,theta):
    x=r*np.sin(theta)*np.cos(phi)
    y=r*np.sin(theta)*np.sin(phi)
    z=r*np.cos(theta)
    return (x, y, z)

def spheric_coordinates_rotate(A,h,A0,h0):
    (x,y,z) = spheric2cart(np.ones(np.size(A)),A, np.pi/2-h)
    dz=np.pi/2-h0;
    
    x_new=x*np.cos(-dz)*np.cos(-A0)-y*np.sin(-A0)*np.cos(-dz)+z*np.sin(-dz);
    y_new=x*np.sin(-A0)+y*np.cos(-A0);
    z_new=-x*np.sin(-dz)*np.cos(-A0)+y*np.sin(-dz)*np.sin(-A0)+z*np.cos(-dz)
    
    r_new, A_new, zen_new =cart2spheric(x_new,y_new,z_new);
    h_new=np.pi/2-zen_new;
    
    return (A_new, h_new)

def hor2pix_keo(Az,Alt,az0,alt0,c,d):
    Azrot, Altrot = spheric_coordinates_rotate(Az,Alt,az0,alt0)
    X=(np.pi/2-Altrot)*np.cos(Azrot)
    Y=(np.pi/2-Altrot)*np.sin(Azrot)
    x=c[0]+c[1]*X+c[2]*Y; # xpix
    y=d[0]+d[1]*X+d[2]*Y; # ypix
    return x, y

def hor2pix_s1c(Az, Alt, az0, alt0, c, d):
    X=-np.cos(Alt)*np.sin(Az-az0)/(np.cos(alt0)*np.cos(Alt)*np.cos(Az-az0)+np.sin(alt0)*np.sin(Alt));
    Y=-(np.sin(alt0)*np.cos(Alt)*np.cos(Az-az0)-np.cos(alt0)*np.sin(Alt))/(np.cos(alt0)*np.cos(Alt)*np.cos(Az-az0)+np.sin(alt0)*np.sin(Alt));
    x=c[0]+c[1]*X+c[2]*Y;
    y=d[0]+d[1]*X+d[2]*Y;
    return x,y

#def bright_relay_keo(flux,n):
#    b=flux*(10**4)/n*3.654288376498914
#    return b
    
#def bright_relay_keo(flux,n):
#    b=flux*(10**4)/n*5.584673336106748
#    return b

def bright_relay_keo(flux,n):
    b=flux*(10**4)/n*3.5719293461623876
    return b

def bright_relay_s1c(flux,n):
    b=flux*(10**4)/n*252.0045968847909
    return b
