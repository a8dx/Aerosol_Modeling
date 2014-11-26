""" 
FILENAME: GISS_MSE_Analysis.py
AUTHOR: Anthony Louis D'Agostino (ald2187@columbia.edu)
LAST EDITED: November 26, 2014
PURPOSE: To manipulate and analyze MSE data from GISS ModelE data retrieved through NASA/GISS FTP server
			to better understand the relationship between MSE and SASM dynamics  
NOTES: GISS Documentation Page: http://data.giss.nasa.gov/modelE/ar5/
DATA: NASA GISS ModelE CMIP5 FTP Repository 

"physics_version=1 (NINT), aerosols and ozone are read in via pre-computed transient aerosol and ozone fields. 
The aerosol indirect effect is parameterized. This corresponds to the rundeck E_AR5_NINT_oR.R" 


########### COMMENTS ###########

# convert climatology function into something more flexible - which includes a time duration argument 
# missing values in Pr series coded as 1e+20 -- need to include checks for this 
# determine how weighting should be included in EOFs and how treated by CDOs 

TASKS
* Take spatial average of each location 
* Specify bounding box outside of the function - make it consistent across all procedures 


* Convert the MSE work into a program with scenario inputs to make it easier to run across experiments 
* Decide whether to keep only .catfiles or also .india - confusion could arise from 
  keeping both?  
* Determine if geographic area for conducting this needs to be changed 
* Make what's common to both bar and scatterplot separate from either, called by both.  


* Review the plot functions I've created and ensure each has a purpose 


VERSION HISTORY:
11/24/2014 - exported MSE and cdo class files to separate files, which are read so long 
             as os directory is consistent 
10/15/2014 - converted MSE calculation and bar plotting from explicit to functions 
9/17/2014 - file created, borrowed largely from pre-existing GISS_Analysis file 


Data Units: 

pr (kg m^-2 s^-1) 
huss (1)
tas (K)

Conversion of given precip units to mm 
  1 kg m-2 s-1 *  60 sec m-1 60 min hr-1 24 hr d-1  = 86,400 x 1 kg m-2 s-1 = 1 mm / day

NASA GISS Model values

Temperature (K)
Specific Humidity (1)

h = (c_p * T) + (L_v * q)
h = [1004 (J/(kg * K)) * T (K) + 2.501e6 * huss (1)]/1000

##
## Example of 2 stream join
##
#  cdo.sub(input = " ".join([ensAvg100, ensAvg108])))  
 


"""

# -- preamble

import datetime as dt 
import glob, os, re  
import calendar 
from cdo import * 
from netCDF4 import Dataset
import pylab as pl 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import ImageGrid
import datetime
from dateutil.relativedelta import *
#from GISS_MSE_Plotting import * 
#from GISS_Analysis_Functions import *    # function helper file 
from cdoCl import cdoClass 
from MSE import mseClass



cdo = Cdo()
cdo.debug = True
EOFswitch = False   # may need to remove this 



def decLinePlot(in_data, titleText, subset, calLength):
  """

  'indata' is a single series in a CDO format, - either from a single pixel
            or spatial average of multiple pixels  
  either from a zonal average or individual location 
  'titleText' is a character string that labels plot title 
  calLength = either 360 or 365 depending on model used 

  """ 
  from itertools import cycle
  import matplotlib.pyplot as plt 

  x,y,z = makeContour(in_data, calLength)    

  # -- total number of decadal periods 
  it = math.floor(z.shape[1]/10)
  indexes = np.array([range(int(it)), range(1,1+int(it))])*10

  max_series = 5   # -- how many series per plot
  num_plot = int(math.floor((indexes.shape[1]-1)/max_series) + 1) 

  # -- initialize output array 
  znew = np.zeros((z.shape[0], int(it)))

  # -- years are extracted from 'y' 
  years = x[0,:]

  # -- vary series markers 
  lines = ["-","--","-.",":"]
  linecycler = cycle(lines)
  decmean = [z[:,indexes[:,i]].mean(axis=1) for i in range(indexes.shape[1])] 
  plt.subplots_adjust(hspace=0.07)

  pnum = 1 
  f = 1 
  x_pr_restrict = [150,275]  # portion kept in final plots 
  x_restrict = [0,calLength]

  # -- take average across individual decades, and label accordingly 
  for c in range(indexes.shape[1]):
    print "Now processing " + str(c) + " of the decadal mean."

    znew[:,c] = decmean[c]

    # -- labeling of correct decades 
    year1 = int(years[0]) + indexes[0,c]
    year2 = int(years[0]) + indexes[1,c]

    print "Year 1: " + str(year1) + " Year 2: " + str(year2)

    ax = plt.subplot(num_plot, 1, pnum)
    box = ax.get_position()

    if pnum == 1:
      ax.set_title(''.join(titleText))

    # -- indexing y is necessary to avoid overplotting, since it's a wide array
    ax.plot(y[:,c], znew[:,c], next(linecycler), label = str(year1)+"-"+str(year2)) 
    #ax.plot(y[:,c], znew[:,c], label = str(year1)+"-"+str(year2))

    if f < max_series and c < range(indexes.shape[1])[-1]: 
      f = f+1
    elif f < max_series and c == range(indexes.shape[1])[-1]:    
      #ax.legend() 

      #box = ax.get_position()
      ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
      # Put a legend to the right of the current axis
      ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

      #ax.legend(shadow=True, fancybox=True)  
      #ax.legend(bbox_to_anchor=(.95, 1.0), loc=2)  
     # plt.ylim(ymin=0)
      ax.set_xlim(x_restrict)
     # ax.set_xlim(y_restrict)
    #  ax.set_ylim([ymin,ymax])
    elif f == max_series:  
      ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])  # Shrink current axis by 20%
      ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox = True) # Put a legend to the right of the current axis
      #ax.legend(shadow=True, fancybox=True)  
      ax.set_xlim(x_restrict)
    #  plt.ylim(ymin=0) 
      plt.tick_params(axis='x', which='both', bottom='off', labelbottom = 'off') 
      #ax1.legend(shadow=True, fancybox=True, loc =2)
    #  ax1.legend(bbox_to_anchor=(1.05, 1), loc=2)
      pnum += 1 
      f = 1   
     # ax.set_xlim(y_restrict)

  #plt.tight_layout()  
  fn = os.path.join(graphics_path, titleText[1] + "_" + titleText[3] + "_" + titleText[5] + "_Decadal_Plot.eps")
  plt.savefig(fn, figsize=(8, 6), dpi = 1200, transparent = True, facecolor='w', edgecolor='k') 
  plt.close()  
  #plt.show()



#ttext = ["Decadal Average ", "MSE ", "for ", titlecase(site), " Under ", str(scenario)]


"""

  fig, (ax1, ax2) = plt.subplots(2, sharex = True, sharey = True)
  fig.suptitle('Contour Plots', fontsize = 14)


  x,y,z = makeContour(pr_var, calLength) 
  im1 = ax1.contourf(x,y,z,ncon, cmap = plt.cm.bwr)
  ax1.set_title('Precipitation')
  divider1 = make_axes_locatable(ax1)
  cax1 = divider1.append_axes("right", size="5%", pad=0.05)
  cbar1 = plt.colorbar(im1, cax=cax1 ,  format="%.3f") # ticks=MultipleLocator(0.2),
  #plt.colorbar(pc, cax=axes)
  

  x,y,z = makeContour(mse_var, calLength)
  im2 = ax2.contourf(x,y,z,ncon, cmap = plt.cm.hot_r)
  ax2.set_title('MSE')
  divider2 = make_axes_locatable(ax2)
  cax2 = divider2.append_axes("right", size = "5%", pad = 0.05)
  cbar2 = plt.colorbar(im2, cax=cax2)  #, ticks=MultipleLocator(0.2), format="%.2f")

  ax2.set_xlabel('Year')

  fig.text(0.04, 0.5, 'Day of Year', va='center', rotation='vertical')
  fig.text(.96, 0.25, 'MSE', va = 'center', rotation = 'vertical')
  fig.text(.96, 0.85, 'Precip', va = 'center', rotation = 'vertical')

  #plt.tight_layout()
  plt.subplots_adjust(top=0.85)
  plt.show()



"""



def makeContour(in_data, calLength):

  """
  Precusor to plotting function - converts a cdo object into a 3d object with z 
  representing the data value.  

  Leap year multiplication is appropriate for GISS output, but may not be the case
  for non-GISS models.  

  calLength: Takes on value of 360 (e.g., HadGEM2) or 365 (e.g., GISS)
  """
  import math 
  from datetime import datetime 

    # -- convert years into integers
  indata = Dataset(in_data)
  yearStr = cdo.showyear(input = in_data)[0].split() 

  # -- create list of strings with date objects
  #true_time = [int(i) for i in cdo.showdate(input = in_data)[0].split()] 
  true_time = cdo.showdate(input = in_data)[0].split()

  # -- take month and year and determine how to cut down to complete calendar years
  month_year_pre = [i[0:7] for i in true_time]
  day_month = [i[5:10] for i in true_time]
  month_year_post = [datetime.strptime(x, '%Y-%m') for x in month_year_pre]

  #jan = [i.month == 1 for i in month_year]
  dec = [i.month == 12 for i in month_year]
  first_jan = [i.month == 1 for i in month_year].index(True)

  if calLength == 365:
    # GISS procedures 
    time_vals = (1461.0/1460.0)*indata.variables['time'][:]  # -- probably unnecessary 
    last_dec = [i for i, j in enumerate(day_month) if j == '12-31'][-1]   
  elif calLength == 360: 
    # HadGEM procedures
    last_dec = [i for i, j in enumerate(day_month) if j == '12-30'][-1]    

  # -- contour data series
  z = indata.variables[cdo.showname(input = in_data)[0]][:][first_jan:last_dec+1]

  # -- can't use strptime for 360 day calendar since Feb 30th gets rejected 
  time_fin = month_year_post[first_jan:last_dec+1]
  yearRange = sorted(list(set([i.year for i in time_fin])))

  zf = z.squeeze().reshape(len(time_fin)/len(yearRange),len(yearRange), order = 'F')

  # -- should auto adjust for 360 or 365 day year 
  time_np = np.array(time_fin)
  dayRange = time_np.reshape(len(yearRange), len(time_fin)/len(yearRange))
  dR = range(dayRange.shape[1])

  X,Y = np.meshgrid(yearRange, dR)
  return X,Y,zf 



"""
makeContour plot repeated if something breaks 


import math  

  # -- convert years into integers
indata = Dataset(in_data)
yearStr = cdo.showyear(input = in_data)[0].split() 

true_time = cdo.showdate(input = in_data)[0].split() 

yearRange = [int(i) for i in yearStr]
time_vals = (1461.0/1460.0)*indata.variables['time'][:]
time_dateFormat = [relativedelta(years = int(yearRange[0]) - 1) + datetime.date.fromordinal(int(round(i))) for i in time_vals]   # convert into date-time format 
# labels = [i.strftime("%B")[:3] + "-" + i.strftime("%d") for i in time_dateFormat]
labels = [i.isoformat() for i in time_dateFormat] 

# -- ensure only full calendar years are included 
n_yr = math.floor(len(time_vals)/ 365)
fullyr_vals = n_yr * 365 

# -- done to ensure the calculations are correct
time_rev = time_vals[0:fullyr_vals]
print "Given number of observations is " + str(len(time_vals)) + " , but concatenating to first " + str(len(time_rev)) + " observations."

# -- contour data series
z = indata.variables[cdo.showname(input = in_data)[0]][:][0:fullyr_vals]
zf = z.squeeze().reshape(len(time_rev)/len(yearRange),len(yearRange), order = 'F')


y = time_rev.reshape(len(yearRange), len(time_rev)/len(yearRange))
dayRange = y[0,:]
X,Y = np.meshgrid(yearRange, dayRange)

day_dateFormat = [relativedelta(years = int(yearRange[0]) - 1) + datetime.date.fromordinal(int(round(i))) for i in dayRange]
#  dayLabels = [i.strftime("%B")[:3] + "-" + i.strftime("%d") for i in day_dateFormat]
dayLabels = [j.isoformat() for j in day_dateFormat] 


return X,Y,zf 

"""



def multiContourPlot(location, pr_var, mse_var, ncon, scen, calLength):
  """
  Purpose: to minimize code reproduction in the creation of contour plots
            using a range of contour band levels.  Does not auto-save plots
            to file, since X-window enables user to save selectively.  

  Variables:           
  location: 
  loctext: string variable used primarily for labeling purposes in output plots 
  var1: first variable to be considered for plotting
  var2: 2nd variable (ideally should create some tuple that spans across
      all variables and automatically generates comparative plots)  
  ncon: number of contour lines in the graph, common to all variables 
  title1/2: Titles correspond to variables being analyzed, original sequence is 
            precip for slot 1 and MSE for slot 2 
  scen: string featuring experiment/scenario (4xCO2, 1%CO2, etc.)           

  """ 

  from mpl_toolkits.axes_grid1 import ImageGrid
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  from matplotlib.ticker import MultipleLocator
  import matplotlib.pyplot as plt 
  from titlecase import titlecase


  #from matplotlib import rc
  #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
  ## for Palatino and other serif fonts use:
  #rc('font',**{'family':'serif','serif':['Palatino']})
  #rc('text', usetex=False)
  #F = plt.figure(1, (5.5, 3.5))

  fig, (ax1, ax2) = plt.subplots(2, sharex = True, sharey = True)
  fig.suptitle('Contour Plots for ' + titlecase(location) + ' Under ' + str(scen), fontsize = 14)

 # plt.rc('text', usetex=True)
 # plt.rc('font', family='serif')

  x,y,z = makeContour(pr_var, calLength) 
  im1 = ax1.contourf(x,y,z,ncon, cmap = plt.cm.bwr)
  ax1.set_title(r'Precipitation (kg m$^2$ s$^{-1}$)')
  divider1 = make_axes_locatable(ax1)
  cax1 = divider1.append_axes("right", size="5%", pad=0.05)
  cbar1 = plt.colorbar(im1, cax=cax1 ,  format="%.3f", ticks=MultipleLocator(0.001))
  #plt.colorbar(pc, cax=axes)
  
  x,y,z = makeContour(mse_var, calLength)
  im2 = ax2.contourf(x,y,z,ncon, cmap = plt.cm.hot_r)
  ax2.set_title(r'MSE (K)')
 # ax2.set_title(r'MSE (  $\circ$ K)')
  divider2 = make_axes_locatable(ax2)
  cax2 = divider2.append_axes("right", size = "5%", pad = 0.05)
  cbar2 = plt.colorbar(im2, cax=cax2, ticks=MultipleLocator(10))  #, format="%.2f")

  ax2.set_xlabel('Year')

  fig.text(0.04, 0.5, r'Day of Year', va='center', rotation='vertical')
  #fig.text(.96, 0.25, 'MSE ' + u"\u00b0" + 'K)', va = 'center', rotation = 'vertical')
  #fig.text(.96, 0.75, 'Precip (kg m$^2$ s$^-1$)', va = 'center', rotation = 'vertical')
  

  fn = os.path.join(graphics_path, scen[3] + "_" + scen[5] + "_MSE_Pr_Contour.eps")
  plt.savefig(fn, figsize=(8, 6), dpi = 1200, transparent = True, facecolor='w', edgecolor='k') 
  #plt.tight_layout()
  plt.subplots_adjust(top=0.85)
  plt.close(fig) 
  #plt.show()



def mse_PlottingAnalysis(path, location, scenario, calLength):

  """
  location: list object containing concluding with "_loc" which has been pre-defined
            with lat and lon values in accordance with cdo requirements for re-mapping 
  scenario: String with experiment specifics (4x, 1%, etc.)          

  ### Trouble with getting "site"_pr to execute, so caved and programmed as following 

  """

  from titlecase import titlecase

  os.chdir(path)

  pr = cdoClass('pr', True)
  pr.analysis() 
  huss = cdoClass('huss', True)
  huss.analysis()
  tas = cdoClass('tas', True)
  tas.analysis() 

  mseTotal = mseClass(tas.post, huss.post) 

  for site in location.keys(): 
    print "Now processing " + str(site)

    sitename = site.split("_")[0]

    # -- get location values for pr and MSE 
    mseTotal.mse_pixel(location[site])
    sitename_mse = mseTotal.pixel  
      
    sitename_pr = cdo.remapnn(location[site], input = pr.post)
    multiContourPlot(sitename, sitename_pr, sitename_mse, 15, ["Decadal Average ", "MSE ", "for ", titlecase(site), " Under ", str(scenario)], 360)
    decLinePlot(sitename_mse, ["Decadal Average ", "MSE ", "for ", titlecase(site), " Under ", str(scenario)], False, calLength)
    decLinePlot(sitename_pr, ["Decadal Average ",  "Precip ", "for ", titlecase(site), " Under ", str(scenario)], False, calLength)





"""
Initial Parameters
"""
#phys313 = 313                     # specify which model physics to work with 
# varSeries = "pr"                  # recall that underscore follows 
#bbox = "58.75,96.25,5,41"         # bbox coordinates applied only if IndiaOnly = True 

#IndiaOnly = True                  # limit grid to aforementioned bbox coordinates
showPlot = False                  # currently only for monthly climatology, but eventually to control all plot outputs
#showEnsembleMembers


"""
Declare path locations 
"""

#basePath = '/Users/xadx/Desktop/Active_Research/Climate_Modelig/GISS_Downloads'
#basePath = '/Volumes/PASSPORTAD/Climate_Modeling/GISS_Downloads/'


modelPath = '/Users/xadx/Desktop/Active_Research/Climate_Modeling'
basePath = os.path.join(modelPath, 'GISS_Downloads') 
graphics_path = os.path.join(modelPath, "Graphics")


# -- GISS models
giss_4x = os.path.join(basePath, "abrupt4xCO2:E2-H_abrupt4xCO2_r1i1p3_day")
giss_rcp85 = os.path.join(basePath, "rcp85:E2-H_rcp85_r2i1p1_day")
giss_e2h_hist = os.path.join(basePath, "E2-H_historical_r6i1p1_day")
giss_e2r_hist  = os.path.join(basePath, "E2-R_historical_r6i1p1_day")

# -- Hadley models
had_base = "/Users/xadx/Desktop/Active_Research/Climate_Modeling/CMIP5_Data/HadGEM2"
had_1pct = os.path.join(had_base, "1pctCO2", "HadGEM2-ES_1pctCO2_r1i1p1_day")
had_4x = os.path.join(had_base, "abrupt4xCO2", "HadGEM2-ES.abrupt4xCO2.day.atmos.day.r1i1p1.v20130214")
rcphad_rcp85 = os.path.join(had_base, "rcp85", "HadGEM2-ES_rcp85_r")    # -- we have 4 model runs for the rcp85 experiment 

# -- have to work on this one!  
# rcp85_dict = {'run1':[rcp85had_rcp]}



# -- coordinates of interest
mumbai_loc = "lon=73.01_lat=19.03"
nagpur_loc = "lon=79.0882_lat=21.1458"
kolkata_loc = "lon=88.3630400_lat=22.5626300" 


# -- location and experiment dictionaries 
location_dict = {'mumbai':mumbai_loc, 
      'nagpur':nagpur_loc,
      'kolkata':kolkata_loc}

hadexp_dict = {'1pct':[had_1pct, "1% CO2 (HadGEM2)"], 
      '4x':[had_4x, "Abrupt 4x CO2 (HadGEM2)"]}

gissexp_dict = {'4x':[giss_4x, "Abrupt 4x CO2 (GISS)"], 
      'RCP85':[giss_rcp85, "RCP 8.5 (GISS)"],
      'e2hhist':[giss_e2h_hist, "Historical (GISS E2-H)"],
      'e2rhist':[giss_e2r_hist, "Historical (GISS E2-R)"]}



"""
Generate MSE values 

mseTotal = mseClass(tas.post, huss.post) 
#  mseTotal.mse_barplot(mumbai_loc, 2908, "Mumbai 2908 4xCO2 Daily MSE")
# mseTotal.mse_barplot(nagpur_loc, 2982, "Nagpur 2982 4xCO2 Daily MSE")

mseTotal.mse_scatterplot(nagpur_loc, (2982,2990), "Nagpur 2982-2990 " + scen + " Daily MSE")

mseTotal.mse_pixel(nagpur_loc)
mseTotal.mse_time_subset(nagpur_loc, (1850, 1900))

# -- would typically just use pixel nameclass, but
# -- because temporal subset, using time_subset instead 
# nagpur_mse = mseTotal.pixel

nagpur_mse = mseTotal.time_subset 

yRange = "1850,1851,1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,1870,1871,1872,1873,1874,1875,1876,1877,1878,1879,1880,1881,1882,1883,1884,1885,1886,1887,1888,1889,1890,1891,1892,1893,1894,1895,1896,1897,1898,1899,1900"
# quick copy and paste job to get it done 
nagpur_pr = cdo.selyear(yRange, input = cdo.remapnn(nagpur_loc, input = pr.post))


mseTotal.mse_pixel(mumbai_loc)
# mseTotal.mse_time_subset(mumbai_loc, (1850, 1900))

mumbai_mse = mseTotal.pixel 
# mumbai_mse = mseTotal.time_subset 
# mumbai_pr = cdo.selyear(yRange, input = cdo.remapnn(mumbai_loc, input = pr.post))

multiContourPlot("nagpur", "pr", "mse", " Precipitation Under "," MSE Under ", scen)
multiContourPlot("mumbai", "pr", "mse", " Precipitation Under "," MSE Under ", scen)

"""


# -- each experiment-location pair for Hadley models
for exp in hadexp_dict.keys():
  mse_PlottingAnalysis(hadexp_dict[exp][0], location_dict, hadexp_dict[exp][1], 360)

# -- each experiment-location pair for GISS models 
for exp in gissexp_dict.keys():
  mse_PlottingAnalysis(gissexp_dict[exp][0], location_dict, gissexp_dict[exp][1], 365)


# -- for each of the RCP8.5 experiments 
# -- ENSEMBLE AVERAGE?? 
for i in range(0,5):
  i_rcp85 = os.path.join(had_rcp85, i, "i1p1")
  mse_PlottingAnalysis(i_rcp85, location_dict, "RCP 8.5 (HadGEM2) Run " + str(i), 360) 





"""
Which day does the peak MSE occur? 
"""


# # -- identify first monsoon day 
# # precip data is in units of 
# pr_convert = 86400
# pr.post = cdo.mulc(pr_convert, input = pr.india)
# pr10 = cdo.gtc('10', input = pr.india, options = "-f nc")  #greater than 10mm 
# pr10 = cdo.mulc('10', input = pr.india, output = ofile, options = "-f nc")  #greater than 10mm 
# monsoon_10mm = cdo.eca_r10mm(input = pr.india)   # identify first date for each pixel that surpasses 10mm threshold 
# m10 = Dataset(monsoon_10mm)
# spacePlot(m10, 0, 'heavy_precipitation_days_index_per_time_period')
