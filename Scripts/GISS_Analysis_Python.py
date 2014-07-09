# GISS_Analysis_Python.py
# Created by: Anthony Louis D'Agostino (ald2187@columbia.edu)
# Last edited: July 9, 2014
# Objective: To manipulate and analyze GISS ModelE data retrieved through NASA/GISS FTP server 
#       later on to develop a script which automates the download of the GISS data to ensure all files are current 


import glob, os, re  
from cdo import * 
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

cdo = Cdo()



cdo.debug = False


def common_prefix(strings):
    """ Find the longest string that is a prefix of all the strings.
    """
    if not strings:
        return ''
    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    return prefix


def eof_standard_plot(obj, value, var, titleText):
  """ Produces some standard plot output and therefore easier to change 
      these settings than manual work. 
      Figure out some way to default to total range for index (i.e., not just 1)
      obj : post-processed netCDF4 object (i.e., not raw original object)
      var : prescribed abbreviation of climate variable of interest 
      value : when working with array, take a single value (i.e., 1st EOF)
  """
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  cliva = obj.variables[var][value,:,:]   # climate variable c
  cliva_units = obj.variables[var].units

  obj.close()  # closes the file

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
  cs = m.pcolor(xi,yi,np.squeeze(cliva))

# Add Grid Lines
  m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
  m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

# Add Coastlines, States, and Country Boundaries
  m.drawcoastlines()
  m.drawstates()
  m.drawcountries()

# Add Colorbar
  cbar = m.colorbar(cs, location='bottom', pad="10%")
  cbar.set_label(cliva_units)

# Add Title
  plt.title(titleText)
  plt.show()      





########### COMMENTS ###########

# not doing a subset anymore 
# only select nc files that have a numeric ending (ensure that generated files do not infringe on this)
# sellonlatbox has explicit bounding region - need to figure out how to introduce variable into this 
# figure out how best to handle the units for precipitation 
# convert climatology function into something more flexible - which includes a time duration argument 
# climatology to be constructed before or after averaging acoss the bounding box?  
# avgArray - need to configure so that can run on more than just a single element within a vector 
# right now just fails if files are not in their place - need to deveop a better file-checking mechanism 
# need to work on ensembleAvg function ASAP 
# missing values in Pr series coded as 1e+20 -- need to include checks for this 
# want to ensure that any array operations are performed AFTER the geographic subset is complete - this minimizes processing requirements 
# currently generating non-melted output in meltObj

# problems registering date range when shifting from E2R to E2H -- needs to be revised in how time is imported 

# Included a time subset into pre-Spatial object to see if EOF actually works 

####### TO-DO ###########


  # generate ensemble mean for comparison purposes 
  # include other climate variables into ST object 
  # remove subset from melt function FILE LIST which restricts number of years analyzed  
  # if stObj works, then "fix" version can be removed 
  # demean by monthly climatology, not total 
  # Need to identify how to better show raster-level data - download raster package.  but requires a different class 
  # could create "India" option which auto-subset to a given bounding box (or another other location)
  # Delete file.list object once function is fully operational 


basePath = '/Users/xadx/Desktop/Active_Research/Climate_Modeling/GISS_Downloads'
EAllDirs = [x[0] for x in os.walk(basePath)]
EHDirs = [f for f in EAllDirs if 'E2-H' in f]



# -- import cell area values 

cellareaRaw = Dataset("areacella_fx_GISS-E2-H-CC_esmFdbk2_r0i0p0.nc", 'r')




# -- see all the paths in vertical list 
print "Only following paths considered in the analysis"
for p in EHDirs:
	print p 

# -- specify which model physics to work with 
modelPhys = 313 
varSeries = "pr"  #recall that underscore follows 
bbox = (58.75,96.25,5,41)   # only if some India flag is turned on 

# -- loops across model output only from a single physics version 
physDirs = [h for h in EHDirs if 'p'+str(modelPhys) in h]
for p in physDirs:
  print "Now processing the following path " + str(p) 
  os.chdir(p) 
  subfiles = glob.glob( varSeries + '_*.nc' )
  filter_function = lambda name: not name.endswith("_cat.nc")  # exlude previously concatenated files 
  goodFiles = filter(filter_function, subfiles)
  print str(len(goodFiles)) + ' files in the path: ' + str(goodFiles)
  
  catPrefix = common_prefix(goodFiles)
  catFileName = catPrefix+'cat'
  print str(goodFiles)

 # if not os.path.exists(os.path.join(p,"outputFiles")):
 #     os.makedirs(os.path.join(p,"outputFiles"))  

  totCat = cdo.cat(input = goodFiles) #,output = os.path.join(p, "output", catFileName + '.nc')    # creates the concatenated file for chosen prefix and modelPhys
  subCat = cdo.sellonlatbox("58.75,96.25,5,41", input = totCat)    #,   input = catFileName + '.nc')
  #subCat = totCat   # leaves room to convert back to the subset version once consistency is achieved 

  cdo.eofspatial(5, input = subCat, output = " ".join(("eof1", "eof2")))   # first 5 spatial EOFs

  eof2 = Dataset('eof2', mode = 'r')

eof_standard_plot(eof2, 1, varSeries, 'EOF1 for ModelPhys = 313')
 







# 	cdo.eoftime(5, input = "/Users/xadx/Desktop/Active_Research/Climate_Modeling/GISS_Downloads/historical:E2-H_historicalMisc_r1i1p105/pr_Amon_GISS-E2-H_historicalMisc_r1i1p105_185001-190012.nc", output = " ".join(("eof1.nc", "eof2.nc")))
# 	cdo.eof(5, input = "/Users/xadx/Desktop/Active_Research/Climate_Modeling/GISS_Downloads/historical:E2-H_historicalMisc_r1i1p105/pr_Amon_GISS-E2-H_historicalMisc_r1i1p105_185001-190012.nc", output = " ".join(("eof1.nc", "eof2.nc")))



#cdo sellonlatbox,58.75,96.25,5,41 output.nc output_subset.nc  # can turn this into a temporary file 


