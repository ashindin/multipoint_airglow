from sympy import *

# init_printing()

fun_str='a1*exp(-(x/a2)**2 - (y/a3)**2)/a4/beta(a5,a6)*((z)/a4)**(a5-1)*((a4-(z))/a4)**(a6-1)'
model_fun_fname='model_drop_fun.py'

fun=sympify(fun_str)

heavy_fun_str1='-z'
heavy_fun_str2='-z+a4'
heavy_fun1=sympify(heavy_fun_str1)
heavy_fun2=sympify(heavy_fun_str2)

x,y,z = symbols('x y z');
xi,eta,zeta = symbols('xi eta zeta')
X, Y, Z = symbols('X Y Z')
a, b, e = symbols('a b e')
R, E, A = symbols('R E A')
phi, lam, h = symbols('phi lam h');
phi_m, lam_m, h_m = symbols('phi_m lam_m h_m');
phi_k, lam_k, h_k = symbols('phi_k lam_k h_k');
a1, a2, a3, a4, a5, a6, a7, a8, a9 = symbols('a1 a2 a3 a4 a5 a6 a7 a8 a9');

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

heavyREA1=heavy_fun1.subs(z,z_fun)
heavyREA2=heavy_fun2.subs(z,z_fun)

heavyREA1_str=str(heavyREA1)
heavyREA2_str=str(heavyREA2)

#heavyREA_str=str(heavyREA).replace('sign','np.signbit')
#print(heavyREA_str)

#fid=open(model_fname,'w')
#fid.write('# sphere model function to integrate\n')
#fid.write(str(funREA)+'\n')
#fid.close()

funREA_str=str(funREA).replace('beta(a5, a6)','beta_a5a6')
#funREA_numpy=funREA_str.replace('exp','np.exp').replace('sin','np.sin').replace('cos','np.cos').replace('sqrt','np.sqrt')

fid=open(model_fun_fname,'w')
fid.write('import numpy as np\n')
fid.write('import numexpr as ne\n')
fid.write('import scipy.special as sp\n')
fid.write('def drop_fun(R,E,A,aa,mod_pos,cam_pos):\n')
fid.write('    a1=aa[0]\n    a2=aa[1]\n    a3=aa[2]\n    a4=aa[3]\n')

fid.write('    if aa[4]<0:\n')
fid.write('        a5=2.\n')
fid.write('        a6=6.\n')
fid.write('    if aa[4]>1.:\n')
fid.write('        a5=6.\n')
fid.write('        a6=2.\n')
fid.write('    if aa[4]<0.5:\n')
fid.write('        a5=2.+aa[4]*8\n')
fid.write('        a6=6.\n')
fid.write('    else:\n')
fid.write('        a5=6.\n')
fid.write('        a6=6.-(aa[4]*8-4.)\n')

fid.write('    beta_a5a6=sp.beta(a5,a6)\n')
fid.write('    norm_coef=((a5 - 1)/(a5 + a6 - 2))**a5*((a6 - 1)/(a5 + a6 - 2))**a6*(a5**2 + 2*a5*a6 - 4*a5 + a6**2 - 4*a6 + 4)/(a4*(a5*a6 - a5 - a6 + 1)*sp.beta(a5, a6))\n')
fid.write('    norm_coef=1/norm_coef\n')
fid.write('    phi_m=mod_pos[0]\n    lam_m=mod_pos[1]\n    h_m=mod_pos[2]\n')
fid.write('    phi_k=cam_pos[0]\n    lam_k=cam_pos[1]\n    h_k=cam_pos[2]\n')
fid.write('    a=6378137.0\n    b=6356752.314245\n    e=0.081819190842965558\n')

fid.write('    NE=ne.evaluate("norm_coef*5.627/10**7*(' + funREA_str + ')")\n')

fid.write('    H1=ne.evaluate("'+ heavyREA1_str +'")\n')
fid.write('    H2=ne.evaluate("'+ heavyREA2_str +'")\n')
fid.write('    H1=np.signbit(H1)\n')
fid.write('    H2=np.signbit(H2)\n')

#fid.write('    return NE*H1*H2\n')
fid.write('    return H1*(~H2)*NE\n')

fid.close()
