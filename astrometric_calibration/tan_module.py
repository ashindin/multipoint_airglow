import numpy as np

def tan_calc_pix2st_coefs(arg, AZ, ALT, x, y):
    
    az0=arg[0]
    alt0=arg[1]
    
    on=np.ones(len(x))
    
    X=-(180/np.pi)*np.cos(ALT)*np.sin(AZ-az0)/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    Y=(180/np.pi)*(np.sin(alt0)*np.cos(ALT)*np.cos(AZ-az0)-np.cos(alt0)*np.sin(ALT))/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    
    matA=np.vstack((on,X,Y)).T
    matB=np.vstack((on,x,y)).T
    
    imatBBB=np.dot(np.linalg.inv(np.dot(matB.T,matB)),matB.T)
    imatAAA=np.dot(np.linalg.inv(np.dot(matA.T,matA)),matA.T)
    
    a=np.dot(imatBBB,X).T
    b=np.dot(imatBBB,Y).T
    c=np.dot(imatAAA,x).T
    d=np.dot(imatAAA,y).T    
    return a,b,c,d
    
def tan_calc_pix2st_coefs_discrep(arg, AZ, ALT, x, y):
    
    az0=arg[0]
    alt0=arg[1]
    
    on=np.ones(len(x))
    
    X=-(180/np.pi)*np.cos(ALT)*np.sin(AZ-az0)/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    Y=(180/np.pi)*(np.sin(alt0)*np.cos(ALT)*np.cos(AZ-az0)-np.cos(alt0)*np.sin(ALT))/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    
    matA=np.vstack((on,X,Y)).T
    matB=np.vstack((on,x,y)).T
    
    imatBBB=np.dot(np.linalg.inv(np.dot(matB.T,matB)),matB.T)
    imatAAA=np.dot(np.linalg.inv(np.dot(matA.T,matA)),matA.T)
    
    a=np.dot(imatBBB,X).T
    b=np.dot(imatBBB,Y).T
    c=np.dot(imatAAA,x).T
    d=np.dot(imatAAA,y).T
    
    x_new=c[0]+c[1]*X+c[2]*Y;
    y_new=d[0]+d[1]*X+d[2]*Y;
    
    sigma=np.sqrt(np.mean((x-x_new)**2+(y-y_new)**2));    
    return sigma
    
def tan_hor2pix(AZ,ALT,az0,alt0,c,d):
    X=-(180/np.pi)*np.cos(ALT)*np.sin(AZ-az0)/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    Y=(180/np.pi)*(np.sin(alt0)*np.cos(ALT)*np.cos(AZ-az0)-np.cos(alt0)*np.sin(ALT))/(np.cos(alt0)*np.cos(ALT)*np.cos(AZ-az0)+np.sin(alt0)*np.sin(ALT));
    X_PIX=c[0]+c[1]*X+c[2]*Y;
    Y_PIX=d[0]+d[1]*X+d[2]*Y;
    return X_PIX, Y_PIX;
    
def tan_pix2hor(x,y,az0,alt0,a,b):
    X=a[0]+a[1]*x+a[2]*y
    Y=b[0]+b[1]*x+b[2]*y
    
    AZ_I=np.arctan2(X,-Y)
    R=np.sqrt(X**2+Y**2)    
    
    ALT_I=np.arctan(180/np.pi/R) # TAN PROJECTION
    
    AZ=az0+np.arctan2(-np.cos(ALT_I)*np.sin(AZ_I),
    np.sin(ALT_I)*np.cos(alt0)-np.cos(ALT_I)*np.sin(alt0)*np.cos(AZ_I))
    
    #if AZ>2*np.pi:
        #AZ-=2*np.pi    
    #if AZ<0:
        #AZ+=2*np.pi
    
    ALT=np.arcsin(np.sin(ALT_I)*np.sin(alt0)+np.cos(ALT_I)*np.cos(alt0)*np.cos(AZ_I))    
    
    return AZ, ALT
