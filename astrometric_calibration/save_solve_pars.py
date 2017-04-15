import numpy as np

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
    a=np.zeros(num_lines)
    b=np.zeros(num_lines)
    c=np.zeros(num_lines)
    d=np.zeros(num_lines)
    for i in range(0,num_lines):
        line_list=lines[i].split(' ')
        fit_names.append(line_list[0])
        az_c[i]=float(line_list[1])
        alt_c[i]=float(line_list[2])
        az0[i]=float(line_list[3])
        alt0[i]=float(line_list[4])
        if len(line_list)==18:
            a[i]=float(line_list[6])
            b[i]=float(line_list[7])
            c[i]=float(line_list[8])
            d[i]=float(line_list[9])
        else:
            a[i]=float(line_list[5])
            b[i]=float(line_list[6])
            c[i]=float(line_list[7])
            d[i]=float(line_list[8])
    return az_c,alt_c,az0,alt0,a,b,c,d
	
def get_med_ind(az_c,alt_c):
    x= np.cos(alt_c)*np.cos(az_c)
    y= np.cos(alt_c)*np.sin(az_c)
    x_med=np.median(x)
    y_med=np.median(y)
    r=np.sqrt((x-x_med)**2+(y-y_med)**2)
#     print(r)
    return np.argmin(r)
