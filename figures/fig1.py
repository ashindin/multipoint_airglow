
# coding: utf-8

# In[1]:

import cartopy.crs as ccrs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# In[7]:

# from cartopy.io.img_tiles import OSM
# import cartopy.feature as cfeature
# import cartopy.io.shapereader as shpreader
from cartopy.io.img_tiles import StamenTerrain


# In[3]:

# from cartopy.feature import NaturalEarthFeature, COLORS


# In[4]:

import cartopy.io.img_tiles as cimgt


# In[5]:

stamen_terrain = StamenTerrain()


# In[13]:

# def scale_bar(ax, length, location=(0.5, 0.05), linewidth=3):
#     """
#     ax is the axes to draw the scalebar on.
#     location is center of the scalebar in axis coordinates ie. 0.5 is the middle of the plot
#     length is the length of the scalebar in km.
#     linewidth is the thickness of the scalebar.
#     """
#     #Projection in metres, need to change this to suit your own figure
#     utm = ccrs.UTM(36)
#     #Get the extent of the plotted area in coordinates in metres
#     x0, x1, y0, y1 = ax.get_extent(utm)
#     #Turn the specified scalebar location into coordinates in metres
#     sbcx, sbcy = x0 + (x1 - x0) * location[0], y0 + (y1 - y0) * location[1]
#     #Generate the x coordinate for the ends of the scalebar
#     bar_xs = [sbcx - length * 500, sbcx + length * 500]
#     #Plot the scalebar
#     ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k', linewidth=linewidth)
#     #Plot the scalebar label
#     ax.text(sbcx, sbcy, str(length) + ' km', transform=utm,
#             horizontalalignment='center', verticalalignment='bottom')


# In[42]:

fig=plt.figure(figsize=(6,3))
# ax = plt.axes(projection=ccrs.Mollweide())
ax = plt.axes(projection=ccrs.PlateCarree())
# ax.coastlines(resolution='10m', zorder=2)
# ax.stock_img()

ax.set_extent((46, 49, 55.5, 56.5))


ax.add_image(stamen_terrain, 10)


plt.plot(46.0991056,56.1434444, marker='o', color='red', markersize=12,
         alpha=0.7)#, transform=ccrs.Geodetic())

plt.text(46.0991056,56.19, u'«Сура» / А',
         verticalalignment='bottom', horizontalalignment='left',
         #transform=text_transform,
         bbox=dict(facecolor='sandybrown', alpha=0.5, boxstyle='round'))

plt.plot(48.7444861,55.9305361, marker='o', color='red', markersize=12,
         alpha=0.7)#, transform=ccrs.Geodetic())

plt.text(48.7444861,55.9305361+0.05, u'Б',
         verticalalignment='bottom', horizontalalignment='left',
         #transform=text_transform,
         bbox=dict(facecolor='sandybrown', alpha=0.5, boxstyle='round'))

ax.gridlines(draw_labels=True)

# scale_bar(ax, 100)

plt.savefig('fig1.eps',dpi=300)
plt.savefig('fig1.pdf',dpi=300)
plt.savefig('fig1.png',dpi=300)

# plt.show()

