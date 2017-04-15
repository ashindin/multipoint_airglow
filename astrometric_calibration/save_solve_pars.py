import numpy as np
import os

def get_solve_coefs(solve_fname):
    fid=open(solve_fname,'r')
    lines=fid.readlines()
    fid.close()
    lines=lines[1::]
    num_lines=len(lines)
    fit_names=[]
    az_c=np.zeros(num_lines)
    alt_c=np.zeros(num_lines)
    az0=np.zeros(num_lines)
    alt0=np.zeros(num_lines)
    a=np.zeros((num_lines,3))
    b=np.zeros((num_lines,3))
    c=np.zeros((num_lines,3))
    d=np.zeros((num_lines,3))
    for i in range(0,num_lines):
        line_list=lines[i].split(' ')
        fit_names.append(line_list[0])
        az_c[i]=float(line_list[1])
        alt_c[i]=float(line_list[2])
        az0[i]=float(line_list[3])
        alt0[i]=float(line_list[4])
        a[i,0]=float(line_list[5])
        a[i,1]=float(line_list[6])
        a[i,2]=float(line_list[7])
        b[i,0]=float(line_list[8])
        b[i,1]=float(line_list[9])
        b[i,2]=float(line_list[10])
        c[i,0]=float(line_list[11])
        c[i,1]=float(line_list[12])
        c[i,2]=float(line_list[13])
        d[i,0]=float(line_list[14])
        d[i,1]=float(line_list[15])
        d[i,2]=float(line_list[16])
    return az_c,alt_c,az0,alt0,a,b,c,d
    
def get_med_ind(az_c,alt_c):
    x= np.cos(alt_c)*np.cos(az_c)
    y= np.cos(alt_c)*np.sin(az_c)
    x_med=np.median(x)
    y_med=np.median(y)
    r=np.sqrt((x-x_med)**2+(y-y_med)**2)
#     print(r)
    return np.argmin(r)

def save_med_solve_pars(solve_fname,save_fname):    
    az_c,alt_c,az0,alt0,a,b,c,d=get_solve_coefs(solve_fname)
    med_ind= get_med_ind(az_c,alt_c)
    az0_med=az0[med_ind]
    alt0_med=alt0[med_ind]
    a_med=a[med_ind,:]
    b_med=b[med_ind,:]
    c_med=c[med_ind,:]
    d_med=d[med_ind,:]
    fid=open(save_fname,'w')
    fid.write("# az0 alt0 a[0] a[1] a[2] b[0] b[1] b[2] c[0] c[1] c[2] d[0] d[1] d[2]\n")
    str_to_file=str(az0_med)+ " " +str(alt0_med)
    str_to_file+=" " + str(a_med[0]) + " " + str(a_med[1]) + " " + str(a_med[2])
    str_to_file+=" " + str(b_med[0]) + " " + str(b_med[1]) + " " + str(b_med[2])
    str_to_file+=" " + str(c_med[0]) + " " + str(c_med[1]) + " " + str(c_med[2])
    str_to_file+=" " + str(d_med[0]) + " " + str(d_med[1]) + " " + str(d_med[2]) + "\n"
    fid.write(str_to_file)
    fid.close()