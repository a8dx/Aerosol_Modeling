import glob, os, re  
from cdo import * 
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

cdo = Cdo()


path = '/Users/xadx/Desktop/Active_Research/Climate_Modeling/GISS_Downloads/historical:E2-H_historicalMisc_r3i1p313'

eof2= Dataset('eof2.nc', mode = 'r')


for v in eof2.variables:
    print v 
    print eof2.variables[v][:]

print eof2.variables['pr'].dimensions


lons = eof2.variables['lon'][:]
lats = eof2.variables['lat'][:]
pr = eof2.variables['pr']
pr_units = eof2.variables['pr'].units

eof2.close()  # closes the file

lon_0 = lons.mean()
lat_0 = lats.mean()

m = Basemap(width=5000000,height=3500000,
            resolution='l',projection='stere',\
            lat_ts=40,lat_0=lat_0,lon_0=lon_0)

# Because our lon and lat variables are 1D, 
# use meshgrid to create 2D arrays 
# Not necessary if coordinates are already in 2D arrays.
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

# Plot Data
cs = m.pcolor(xi,yi,np.squeeze(pr))

# Add Grid Lines
#m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
#m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
m.drawcoastlines()
m.drawstates()
m.drawcountries()

# Add Colorbar
cbar = m.colorbar(cs, location='bottom', pad="10%")
cbar.set_label(pr_units)

# Add Title
plt.title('First EOF')

plt.show()
