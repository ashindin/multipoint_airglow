
# coding: utf-8

# In[1]:

from sympy import *
# init_printing()


# In[2]:

filename='coord_sys_replacement.dat'


# In[3]:

x,y,z = symbols('x y z');
xi,eta,zeta = symbols('xi eta zeta') 
X, Y, Z = symbols('X Y Z')
a, b, e = symbols('a b e')
R, E, A = symbols('R E A')
phi, lam, h = symbols('phi lam h');
phi_m, lam_m, h_m = symbols('phi_m lam_m h_m');
phi_k, lam_k, h_k = symbols('phi_k lam_k h_k');
a1, a2 = symbols('a1 a2');


# In[4]:

ROT_M=Matrix([[-sin(lam),-cos(lam)*sin(phi),cos(lam)*cos(phi)],
          [cos(lam),-sin(lam)*sin(phi),sin(lam)*cos(phi)],
          [0,cos(phi),sin(phi)]])
# pprint(ROT_M)
# pprint(R.T)


# In[5]:

N_m=a/(sqrt(1-e**2*sin(phi_m)**2))
X_m=(N_m+h_m)*cos(phi_m)*cos(lam_m)
Y_m=(N_m+h_m)*cos(phi_m)*sin(lam_m)
Z_m=(b**2/a**2*N_m+h_m)*sin(phi_m)

N_k=a/(sqrt(1-e**2*sin(phi_k)**2))
X_k=(N_k+h_k)*cos(phi_k)*cos(lam_k)
Y_k=(N_k+h_k)*cos(phi_k)*sin(lam_k)
Z_k=(b**2/a**2*N_k+h_k)*sin(phi_k)
# pprint(Z_k)

dX_km=X_k-X_m
dY_km=Y_k-Y_m
dZ_km=Z_k-Z_m


# In[6]:

ENU_coord=ROT_M.T.subs([(phi,phi_m),(lam,lam_m)])*Matrix([[X],[Y],[Z]])
# pprint(ENU_coord)


# In[7]:

ECEF_coord=ROT_M.subs([(phi,phi_k),(lam,lam_k)])*Matrix([[xi],[eta],[zeta]])
# pprint(ECEF_coord)


# In[8]:

ENU2_coord=ENU_coord.subs([(X,ECEF_coord[0]+dX_km),(Y,ECEF_coord[1]+dY_km),(Z,ECEF_coord[2]+dZ_km)])
ENU2_coord=simplify(ENU2_coord)
# pprint(dZ_km.subs([(h_k,0),(h_m,0),(phi_m,-pi/2),(phi_k,pi/2)])) # should be > 0
# pprint(ENU2_coord)


# In[9]:

ELAZ_coord=ENU2_coord.subs([(xi,R*cos(E)*sin(A)),(eta,R*cos(E)*cos(A)),(zeta,R*sin(E))])


# In[10]:

x_new=collect(ELAZ_coord[0], R)
y_new=collect(ELAZ_coord[1], R)
z_new=collect(ELAZ_coord[2], R)


# In[11]:

fid=open(filename,'w')
fid.write('# [x, y, z] replacement to R, E, A coordinates:\n')
fid.write(str(x_new)+'\n')
fid.write(str(y_new)+'\n')
fid.write(str(z_new)+'\n')
fid.close()

