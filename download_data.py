
import requests
import json
import os

def download_check_unpack_data():

	urls_filename="./data/download_urls.txt"
	fid=open(urls_filename,'r')
	urls=fid.readlines()
	fid.close()

	for i in range(len(urls)):
		if urls[i][-1]=='\n':
			urls[i]=urls[i][:-1]

	def get_direct_url(url_yandex):
		r = requests.get("https://cloud-api.yandex.net/v1/disk/public/resources/download?public_key="+url_yandex)
		if r.status_code==200:
			return json.loads(r.content)['href']
		else:
			return 1
	direct_urls=[]
	for i in range(len(urls)):
		direct_url=get_direct_url(urls[i])
		if direct_url==1:
			print("Can not get direct url for " + urls[i])
			break
		direct_urls.append(direct_url)
	# direct_urls
	filenames=["./data/"+s.split("filename=")[-1].split("&")[0] for s in direct_urls]


	for i in range(len(direct_urls)):
		if os.path.exists(filenames[i])==True:
			print("File " + filenames[i] + " already exists")
			continue
		else:
			print("Downloading file " + filenames[i])
			os.system('wget -O ' + filenames[i] +' "'+ direct_urls[i]+'"')

	if os.system("cd data && md5sum -c twopoint_data_optic_2014_2016.md5")==0:

		directory="./data/140824/keo"
		if not os.path.exists(directory):
			os.makedirs(directory)
		ret=os.system("tar xvjf "+filenames[0] +" -C "+ directory)
		
		directory="./data/140824/s1c"
		if not os.path.exists(directory):
			os.makedirs(directory)
		ret+=os.system("tar xvjf "+filenames[1] +" -C "+ directory)
		
		directory="./data/140824/photometers"
		if not os.path.exists(directory):
			os.makedirs(directory)
		ret+=os.system("tar xvjf "+filenames[2] +" -C "+ directory)
		
		directory="./data/140826/s1c"
		if not os.path.exists(directory):
			os.makedirs(directory)
		ret+=os.system("tar xvjf "+filenames[3] +" -C "+ directory)
		
		directory="./data/140826/photometers"
		if not os.path.exists(directory):
			os.makedirs(directory)
		ret+=os.system("tar xvjf "+filenames[4] +" -C "+ directory)
		
		directory="./data/140826/keo"
		if not os.path.exists(directory):
			os.makedirs(directory)
		ret+=os.system("tar xvjf "+filenames[5] +" -C "+ directory)
	  
		directory="./data/160829/keo"
		if not os.path.exists(directory):
			os.makedirs(directory)
		ret+=os.system("tar xvjf "+filenames[6] +" -C "+ directory)
			
		directory="./data/160829/sbig"
		if not os.path.exists(directory):
			os.makedirs(directory)
		ret+=os.system("tar xvjf "+filenames[7] +" -C "+ directory)
		
		directory="./data/160829/hf"
		if not os.path.exists(directory):
			os.makedirs(directory)
		ret+=os.system("tar xvjf "+filenames[8] +" -C "+ directory)
		
		if ret==0:
			return 0
		else:
			print("Error! Something wrong with data unpacking!")
			return 1
		
	else:
		print("Downloaded files are corrupted!")
		return 1

ret=download_check_unpack_data()
if ret==0:
	print("Data are downloaded, checked and unpacked! All is fine!")