
# coding: utf-8

# In[5]:

import requests
import json
import os


# In[2]:

urls_filename="./data/download_urls.txt"
fid=open(urls_filename,'r')
urls=fid.readlines()
fid.close()
for i in range(len(urls)):
    if urls[i][-1]=='\n':
        urls[i]=urls[i][:-1]


# In[ ]:

# urls


# In[4]:

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

# In[ ]:

for i in range(len(direct_urls)):
    print("Downloading file " + str(i))
    os.system('wget -O ' + filenames[i] +' "'+ direct_urls[i]+'"')

