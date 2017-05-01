import os
import subprocess
import sys

def main_script(python_cmd):
    os.system(python_cmd+" make_master_darks.py")
    os.system(python_cmd+" s1c_sp_calibration.py")
    os.system(python_cmd+" keo_sp_calibration.py")
    os.system(python_cmd+" sbig_sp_calibration.py")

os.chdir('./spectrophotometric_calibration/')

python_cmd="python"

try:
    CP=subprocess.run([python_cmd, "--version"], check=True, stderr=subprocess.PIPE)
except FileNotFoundError:
    python_cmd="python3"
    CP=subprocess.run([python_cmd, "--version"], check=True, stderr=subprocess.PIPE)

python_version_str=str(CP.stderr)[9:14]

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
