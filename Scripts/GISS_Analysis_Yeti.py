""" 
FILENAME: GISS_Analysis_Python.py

AUTHOR: Anthony Louis D'Agostino (ald2187@columbia.edu)

LAST EDITED: July 15, 2014

PURPOSE: To manipulate and analyze GISS ModelE data retrieved through NASA/GISS FTP server 
    later on to develop a script which automates the download of the GISS data to ensure all files are current 


NOTES

GISS Documentation Page: http://data.giss.nasa.gov/modelE/ar5/

"physics_version=1 (NINT), aerosols and ozone are read in via pre-computed transient aerosol and ozone fields. 
The aerosol indirect effect is parameterized. This corresponds to the rundeck E_AR5_NINT_oR.R" 


* Need to complete monthlyPlot function 
* How to restrict analysis only to last 30 years (something along these lines)
* Identify the cause for the month 1 allocation error in ymonsub
Q: Take the rootgrps approach (from netcdf4 docs) to compile all objects into a single groupset? 
Have hard-coded 1850 as the starting year for the simplePlot -- needs to be fixed 
* Can later convert spEOF into more all-purpose, for any type of EOF 
Learn more about how the aerosols are prescribed and information on their concentration levels 
* Figure out how to make the ensemble average inside the genEnsemble func
# shortened the period of time under simplePlot for all 3 variables
* Need to interpret WMGHG from GISS
* Ensure that my interpretation of the various GISS scenarios is accurate 
* need to do anything to differentiate between historical and historical2 paths?? 

########### COMMENTS ###########

# only select nc files that have a numeric ending (ensure that generated files do not infringe on this)
# figure out how best to handle the units for precipitation 
# convert climatology function into something more flexible - which includes a time duration argument 
# missing values in Pr series coded as 1e+20 -- need to include checks for this 

TASKS
Make function to read ensemble average 


VERSION HISTORY:
7/15/14 - converting ensemble averaging into a function to be called for individual model physics 
          able to successfully create monthly climatology which needs to be verified 
"""

# -- preamble

import datetime as dt 
import glob, os, re  
import calendar 
from cdo import * 
from netCDF4 import Dataset
from pylab import * 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import ImageGrid
from GISS_Analysis_Functions import *    # function helper file 

cdo = Cdo()

cdo.debug = False
EOFswitch = False   # may need to remove this 


#### FUNCTIONS #####

def ncdump(nc_fid):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.  ******
    From: http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html
    '''
    nc_dims = [dim for dim in nc_fid.dimensions]  # gets your nc dimensions
    nc_vars = [var for var in nc_fid.variables]  # gets your nc variables
    # NetCDF file description
    print "NetCDF file description:\n\t", repr(nc_fid.description)
    # Dimension shape information.
    print "NetCDF dimension information:"
    for dim in nc_dims:
        print '\tName:', dim, "Size:", len(nc_fid.dimensions[dim])
        for ncattr in nc_fid.variables[dim].ncattrs():
            print '\t\t', ncattr,\
            repr(nc_fid.variables[dim].getncattr(ncattr))
    # Variable information.
    print "NetCDF variable information:"
    for var in nc_vars:
        if var not in nc_dims:
            print '\tName:', var, "Size:", nc_fid.variables[var].size,\
                  "Dimensions:", nc_fid.variables[var].dimensions
            for ncattr in nc_fid.variables[var].ncattrs():
                print '\t\t', ncattr,\
                      repr(nc_fid.variables[var].getncattr(ncattr))
    return nc_dims, nc_vars




def eof_standard_plot(object, value, var, titleText, filename):
  """ Produces some standard plot output and therefore easier to change 
      these settings than manual work. 
      Figure out some way to default to total range for index (i.e., not just 1)
      obj : post-processed netCDF4 object (i.e., not raw original object)
      var : prescribed abbreviation of climate variable of interest 
      value : when working with array, take a single value (i.e., 1st EOF)
  """
  obj = Dataset(object, mode = 'r')
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  cliva = obj.variables[var][value,:,:]   # climate variable
  cliva_units = obj.variables[var].units

  

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
  #plt.show()
  #plt.draw()   
  fn = str(filename) + ".eps"     # hard-coding here can be improved
  plt.savefig(fn)
  #plt.show()
  #fig.show() 
  #def ensemble_avg(modphys)
  obj.close()  # closes the file



def simplePlot(obj):
  """
    Highly convertible definition to produce plots on the fly from netCDF4 objects 
    and user-specified climate variables. 

    Assumes that yvar is in a (N,) array format (i.e., 1D object) 
    obj:  post-'Dataset' object (until figuring out how to Dataset inside function)
    time: variable plotted on x-axis
    yvar: variable plotted on y-axis
    lab:  plot label text  
  """

  yvar = "pr"
  lab = "Insert text here"
  time = obj.variables['time'][:]
  dt_time = [dt.date(1850,1,1) + dt.timedelta(days=t) for t in time]
  #time.shape = time.shape[0]    # keep first dimension only
  yt = obj.variables[yvar][:]
  yt.shape = yt.shape[0]    # keep first dimension only
  yt_units = obj.variables[yvar].units
  #ticks = num2date(t1, "days since 1850-01-01 00:00:00")
  #fig = figure()
  plot(dt_time,yt, label = lab)
  xlabel("Time")
  xticks(rotation = 55)
  #autofmt_xdate()
  show()
  obj.close()


def genEnsemble(modelphys, plotMembersClima = False):
  """
  Given a model physics numeric value, provides the ensemble average with an option to plot 
  a comparative climatology for all ensemble members.  Output is not demeaned, enabling easier 
  comparison with output from other model physics runs.  
  """

  # -- see all the paths in vertical list 
  EAllDirs = [x[0] for x in os.walk(basePath)]
  EHDirs = [f for f in EAllDirs if 'E2-H' in f]
  print "Only following paths considered in the analysis"
  for p in EHDirs:
    print p 

  count = 0   # not certain what role of this is   
  physDirs = [h for h in EHDirs if 'p'+str(modelphys) in h]   # only paths with specified model physics are considered 
  ensembleTotal = [] 
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

    totalCat = cdo.cat(input = goodFiles) #, options = '-f nc') #,output = os.path.join(p, "output", catFileName + '.nc')    # creates the concatenated file for chosen prefix and modelPhys
    
    if (IndiaOnly == True):
      print "Restricting coordinates to " + str(bbox)
      totalCat = cdo.sellonlatbox(str(bbox), input = totalCat) #, options =  '-f nc')    #,   input = catFileName + '.nc')
    else: 
      print "Maintaining total coordinates!"

    ensembleTotal.append(totalCat)
    count += 1    # identifier for ensemble total 

    if (plotMembersClima == True):
      ensemblePlot()

    print str(len(ensembleTotal)) + " length of ensemble total"  



    # alt include ensAverage = cdo...
  print ensembleTotal  
  #ensA = cdo.ensavg(input = ensembleTotal)  #, options =  '-f nc')         # ensemble average 
  #return ensAverage 
  return ensembleTotal





def monthlyPlot(object, var, title):
  
  # time goes 0:11

  obj = Dataset(object)
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  cliva_units = obj.variables[var].units

  lon_0 = lons.mean()
  lat_0 = lats.mean()
  

  fig, axes = plt.subplots(nrows = 4, ncols = 3)

  count = 0    # to identify slices 
  
  grid = AxesGrid(fig, [0.05,0.01,0.9,0.9], # similar to subplot(132)
                    nrows_ncols = (4, 3),
                    axes_pad = 0.25,
                    cbar_mode='single',
                    label_mode = "L",
                    cbar_location = "bottom",
                    share_all=True
                    )


  for ax in axes.flat:
    
    cliva = obj.variables[var][count,:,:]   # climate variable

    #m_ax = Basemap(ax = ax, width=5000000,height=3500000,resolution='l',projection='stere', lat_ts=40,lat_0=lat_0,lon_0=lon_0)
    m_ax = Basemap(ax = ax, width=5000000,height=3500000, projection = 'ortho', lat_ts = 40, lat_0 = lat_0, lon_0 = lon_0)

    # Because our lon and lat variables are 1D, 
    # use meshgrid to create 2D arrays 
    # Not necessary if coordinates are already in 2D arrays.
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m_ax(lon, lat)

    cs = m_ax.pcolor(xi,yi,np.squeeze(cliva))


    #cbar = m_ax.colorbar(cs, location='bottom', pad="40%")

    m_ax.drawcoastlines()
    m_ax.drawcountries()
    
      # Add Grid Lines
    m_ax.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m_ax.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)
    #ax.set_title('title')

    plt.title(calendar.month_abbr[count+1])
    count += 1 
    #plt.show()
    #cbar.set_label(cliva_units)


  cbar = fig.colorbar(cs, cax = grid.cbar_axes[0], orientation='horizontal')  
  plt.suptitle(title)
  #plt.tight_layout()
  #plt.subplots_adjust(left = 0.1, right = 0.1, top = 0.1, bottom = 0.1)
  plt.show()   



def spacePlot(obj, slice, var):
  """
  Will assume post-Dataset for now, until can get working 
  """
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  cliva = obj.variables[var][value,:,:]   # climate variable
  cliva_units = obj.variables[var].units

  

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
  #plt.show()
  #plt.draw()   
  fn = "EOF1" + str(goodFiles[0][0:-3]) + ".eps"     # hard-coding here can be improved
  plt.savefig(fn)
  #plt.show()
  #fig.show() 
  #def ensemble_avg(modphys)
  obj.close()  # closes the file  



def ensemblePlot():
  """
  ensemblePlot: Appends monthly climatologies for each ensemble member 
  as looping through a vector with concatenated member files. 
  Currently requires no input since it feeds off output objects from 
  standard `genEnsembleAvg' function.
  """ 
  ensNumber = 1    # extract ensemble member number from the file path eventually 
  for member in ensembleTotal:
    y = cdo.ymonavg(input = "-fldavg " + member, options = '-f nc')
    df = Dataset(y, mode = 'r')  
    dfTime = df.variables['time'][:]    #climatology time 
    dfVar = df.variables[varSeries][:]
    dfVar.shape = (12)   # reshape to 1D 
    seriesLab = str(modelPhys) + "Member " + str(ensNumber)
    plot(dfT, dfVar, label = seriesLab)
    ensNumber += 1 
    legend(loc = 'upper left')
    show() 


def spEOF(obj, number):
  """
  spEOF:  
    Returns first 'number' of spatial EOFs from obj.  
    obj: netCDF (post-Dataset?) object from which to create EOFs
    number: first N number EOFs are returned
  """
  cdo.eofspatial(number, input = obj, output = " ".join(("eof1", "eof2")))   # first 5 spatial EOFs
  eof2 = Dataset('eof2', mode = 'r')
  eof_standard_plot(eof2, 1, varSeries, 'EOF1 for ModelPhys = 313')



# -- Initial Parameters
phys313 = 313                     # specify which model physics to work with 
varSeries = "pr"                  # recall that underscore follows 
bbox = "58.75,96.25,5,41"         # bbox coordinates applied only if IndiaOnly = True 
IndiaOnly = True                  # limit grid to aforementioned bbox coordinates
showPlot = False                  # currently only for monthly climatology, but eventually to control all plot outputs
#showEnsembleMembers


#basePath = '/Users/xadx/Desktop/Active_Research/Climate_Modeling/GISS_Downloads'
basePath = '/vega/sscc/work/users/ald2187/GISS_Downloads'

# -- import cell area values 

#cellareaRaw = Dataset("areacella_fx_GISS-E2-H-CC_esmFdbk2_r0i0p0.nc", 'r')



# -- 313 physics    : 13  Long-lived GHGs (CO2, N2O, small miscellaneous gases, but not CH4) 
#ensAvg313 = cdo.ensavg(input = genEnsemble(313))
#moCli313 = cdo.ymonavg(input = ensAvg313)

# -- 300 physics    : 00  all forcings except aerosol indirect effects (pN00 is equivalent to pN in this experiment)




"""
Configuration of models under various physics.  
"""


# # -- 100 physics     : 00  all forcings except aerosol indirect effects (pN00 is equivalent to pN in this experiment)
# ensAvg100 = cdo.ensavg(input = genEnsemble(100))
# moCli100 = cdo.ymonavg(input = ensAvg100)
# deMeaned100 = cdo.ymonsub(input = " ".join([ensAvg100, cdo.ymonavg(input = ensAvg100)]), options = '-f nc')    # demeaning monthly climatology by pixel
# print "\n Successfully processed models under physics 100. \n  all forcings except aerosol indirect effects (pN00 is equivalent to pN in this experiment) \n"


# # -- 106 physics    : 06  anthro tropospheric aerosol (direct effect) only (conc)
# ensAvg106 = cdo.ensavg(input = genEnsemble(106))
# moCli106 = cdo.ymonavg(input = ensAvg106)
# deMeaned106 = cdo.ymonsub(input = " ".join([ensAvg106, cdo.ymonavg(input = ensAvg106)]), options = '-f nc')    # demeaning monthly climatology by pixel
# print "\n Successfully processed models under physics 106. \n anthro tropospheric aerosol (direct effect) only (conc) \n "

# diff106 = cdo.sub(input = " ".join([ensAvg100, ensAvg106]))
# cliDiff106 = cdo.ymonavg(input = diff106)
# # summer months only 


# spDiffCli106 = cdo.ymonavg(input = "-fldavg " + cdo.ymonavg(input = diff106))   # spatial average of 
# x = Dataset(spDiffCli106)
# simplePlot(x)



# # -- 107 physics     : 07   anthro tropospheric aerosol (direct and indirect effects) only (conc)
# ensAvg107 = cdo.ensavg(input = genEnsemble(107))
# moCli107 = cdo.ymonavg(input = ensAvg107)
# deMeaned107 = cdo.ymonsub(input = " ".join([ensAvg107, cdo.ymonavg(input = ensAvg107)]), options = '-f nc')    # demeaning monthly climatology by pixel
# print "\n Successfully processed models under physics 107. \n anthro tropospheric aerosol (direct and indirect effects) only (conc) \n"


# -- 108 physics     : 08  BC on snow only
# ensAvg108 = cdo.ensavg(input = genEnsemble(108))
# moCli108 = cdo.ymonavg(input = ensAvg108)
# deMeaned108 = cdo.ymonsub(input = " ".join([ensAvg108, cdo.ymonavg(input = ensAvg108)]), options = '-f nc')    # demeaning monthly climatology by pixel
# print "\n Successfully processed models under physics 108. \n  BC on snow only.  \n "


# -- 109 physics     
ensAvg109 = cdo.ensavg(input = genEnsemble(109))
moCli109 = cdo.ymonavg(input = ensAvg109)
deMeaned109 = cdo.ymonsub(input = " ".join([ensAvg109, cdo.ymonavg(input = ensAvg109)]), options = '-f nc')    # demeaning monthly climatology by pixel
print "\n Successfully processed models under physics 109. \n  SPECIFY HERE  \n "

eof_standard_plot(ensAvg109, 1, "pr", "Model Physics 109 Ensemble: Precipitation EOF1", "mp309EOF1")

# -- 110 physics     : 10  anthro tropospheric aerosol (via emissions of SO2, BC, OC, NH3) 





# -- Aerosol only effects 

# ensure common time series 












# double check to make sure differencing is in right direction   

# -- manually determine if this is conducting the operation properly 
#monSp = cdo.fldmean(input = deMeaned)
#monS = Dataset(monSp)
#simplePlot(monS)



#first convert from netcdf temp file to run simplePlot, o/w "global name 'Dataset' is not defined"

#x = Dataset(spAvgDemeaned, mode = 'r')
#simplePlot(x)


#files = [ensembleTotal, ensAverage]
#g = cdo.sub(input=" ".join((ensembleTotal, ensAverage)), options = '-f nc')

print "Script successfully completed!"
#plt.show()  
#plt.savefig("fake_eof.png")

