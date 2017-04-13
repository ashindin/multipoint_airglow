filename = "cygwin_install_astrometry.sh"
fileContents = open(filename,"r").read()
f = open(filename,"w", newline="\n")
f.write(fileContents)
f.close()

filename = "./data/twopoint_data_optic_2014_2016.md5"
fileContents = open(filename,"r").read()
f = open(filename,"w", newline="\n")
f.write(fileContents)
f.close()

filename = "./data/download_urls.txt"
fileContents = open(filename,"r").read()
f = open(filename,"w", newline="\n")
f.write(fileContents)
f.close()
