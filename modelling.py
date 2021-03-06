import os
import subprocess
import sys

def main_script(python_cmd):
    os.system(python_cmd+" coord_sys_gen.py")
    os.system(python_cmd+" sphere_model_gen.py")
    os.system(python_cmd+" spheroid_model_gen.py")
    os.system(python_cmd+" ellipsoid_model_gen.py")
    os.system(python_cmd+" spheroid_incl_model_gen.py")
    os.system(python_cmd+" drop_model_gen.py")

os.chdir('./modelling')

python_cmd="python"

try:
    CP=subprocess.run([python_cmd, "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except FileNotFoundError:
    python_cmd="python3"
    CP=subprocess.run([python_cmd, "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

python_version_str=str(CP.stdout + CP.stderr)[9:14]

main_v=int(python_version_str.split('.')[0])
major_v=int(python_version_str.split('.')[1])
minor_v=int(python_version_str.split('.')[2])

if main_v==3 and major_v>=6:
    print("You have a correct python version (>=3.6)")
    main_script(python_cmd)
    sys.exit()
else:
    python_cmd="python3"
    try:
        CP=subprocess.run([python_cmd, "--version"], check=True, stderr=subprocess.PIPE)
    except FileNotFoundError:
        print("You don't have a correct python version")

        if main_v==3 and major_v>=6:
            print("You have a correct python version")
            main_script(python_cmd)
            sys.exit()
        else:
            print("You don't have a correct python version (>=3.6)")
