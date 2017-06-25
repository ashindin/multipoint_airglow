from sympy import *
# init_printing()

fun_str='a1*exp(-(x**2 + y**2 + z**2) /a2**2)'
model_fname='sphere.dat'
model_fun_fname='model_sphere_fun.py'

fun=sympify(fun_str)

x,y,z = symbols('x y z');
xi,eta,zeta = symbols('xi eta zeta')
X, Y, Z = symbols('X Y Z')
a, b, e = symbols('a b e')
R, E, A = symbols('R E A')
phi, lam, h = symbols('phi lam h');
phi_m, lam_m, h_m = symbols('phi_m lam_m h_m');
phi_k, lam_k, h_k = symbols('phi_k lam_k h_k');
a1, a2 = symbols('a1 a2');

coord_replacement_fn='coord_sys_replacement.dat'


fid=open(coord_replacement_fn,'r')
fid.readline()
x_str=fid.readline()[0:-1]
y_str=fid.readline()[0:-1]
z_str=fid.readline()[0:-1]
fid.close()

x_fun=sympify(x_str)
y_fun=sympify(y_str)
z_fun=sympify(z_str)

funREA=fun.subs([(x,x_fun),
                (y,y_fun),
                (z,z_fun)])

fid=open(model_fname,'w')
fid.write('# sphere model function to integrate\n')
fid.write(str(funREA)+'\n')
fid.close()

funREA_str=str(funREA)
funREA_numpy=funREA_str.replace('exp','np.exp').replace('sin','np.sin').replace('cos','np.cos').replace('sqrt','np.sqrt')

fid=open(model_fun_fname,'w')
fid.write('import numpy as np\n')
fid.write('def sphere_fun(R,E,A,aa,mod_pos,cam_pos):\n')
fid.write('    a1=aa[0]\n    a2=aa[1]\n')
fid.write('    phi_m=mod_pos[0]\n    lam_m=mod_pos[1]\n    h_m=mod_pos[2]\n')
fid.write('    phi_k=cam_pos[0]\n    lam_k=cam_pos[1]\n    h_k=cam_pos[2]\n')
fid.write('    a=6378137.0\n    b=6356752.314245\n    e=0.081819190842965558\n')
fid.write('    return ' + funREA_numpy + '\n')
fid.close()
