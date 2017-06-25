from sympy import *
# init_printing()

fun_str='a1*exp(-(x**2 + y**2 + z**2) /a2**2)'
model_fname='sphere.dat'

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
