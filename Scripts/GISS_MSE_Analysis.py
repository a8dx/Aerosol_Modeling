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

# modified pixel -> spatial average in both MSE (self.pixel) and in the mse plotting analysis function 
# convert climatology function into something more flexible - which includes a time duration argument 
# missing values in Pr series coded as 1e+20 -- need to include checks for this 
# determine how weighting should be included in EOFs and how treated by CDOs 

TASKS
* Work on the makeBoxPlot function -- 
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

import glob, os, re  
import calendar 
from cdo import * 
from netCDF4 import Dataset
import pylab as pl 
import numpy as np
from itertools import cycle
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from datetime import datetime 
from dateutil.relativedelta import *
#from GISS_MSE_Plotting import * 
#from GISS_Analysis_Functions import *    # function helper file 
from cdoCl import cdoClass 
from MSE import mseClass
import math 
import pandas as pd
from titlecase import titlecase




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

  # -- need range to get full interval, not just end points
  decmean = [z[:,range(indexes[:,i][0],indexes[:,i][1])].mean(axis=1) for i in range(indexes.shape[1])] 
  plt.subplots_adjust(hspace=0.07)

  pnum = 1 
  f = 1 
  x_pr_restrict = [150,275]  # portion kept in final plots 
  x_restrict = [0,calLength]

  # -- take average across individual decades, and label accordingly 
  for c in range(indexes.shape[1]):
    print "Now processing " + str(c) + " of the decadal mean."

    znew[:,c] = decmean[c]

    # -- labeling of correct decades, intervals are decades  
    year1 = int(years[0]) + indexes[0,c] -1
    year2 = int(years[0]) + indexes[1,c] -1  

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
  labels = [item.get_text() for item in ax.get_xticklabels()]

  fn = os.path.join(graphics_path, titleText[5], titleText[1] + "_" + titleText[3] + "_" + titleText[5] + "_Decadal_Plot.eps")
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
  dec = [i.month == 12 for i in month_year_post]
  first_jan = [i.month == 1 for i in month_year_post].index(True)

  if calLength == 365:
    # GISS procedures 
    time_vals = (1461.0/1460.0)*indata.variables['time'][:]  # -- probably unnecessary 
    last_dec = [i for i, j in enumerate(day_month) if j == '12-31'][-1]   
  elif calLength == 360: 
    # HadGEM procedures (30 day months)
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




  #from matplotlib import rc
  #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
  ## for Palatino and other serif fonts use:
  #rc('font',**{'family':'serif','serif':['Palatino']})
  #rc('text', usetex=False)
  #F = plt.figure(1, (5.5, 3.5))

  fig, (ax1, ax2) = plt.subplots(2, sharex = True, sharey = True)
  fig.suptitle('Contour Plots for ' + titlecase(location) + ' Under ' + str(scen[4]), fontsize = 14)

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
  

  fn = os.path.join(graphics_path, scen[5], scen[3] + "_" + scen[5] + "_MSE_Pr_Contour.eps")
  plt.savefig(fn, figsize=(8, 6), dpi = 1200, transparent = True, facecolor='w', edgecolor='k') 
  #plt.tight_layout()
  plt.subplots_adjust(top=0.85)
  plt.close(fig) 
  #plt.show()



def makeMSE(path, location, scenario, calLength):

  """
  Broken out to access MSE objects without going through plotting loops 
  10/30 used to determine how best to merge muiltiple pixels into a single analysis 

  does remapnn do the same thing as taking a field average across space? 

  location: list object containing concluding with "_loc" which has been pre-defined
            with lat and lon values in accordance with cdo requirements for re-mapping 
  scenario: String with experiment specifics (4x, 1%, etc.)    

  Returns an India-subsetted precip file as well as a (India-subsetted?) MSE file
  """

  os.chdir(path)

  pr = cdoClass('pr', True)
  pr.analysis() 
  huss = cdoClass('huss', True)
  huss.analysis()
  tas = cdoClass('tas', True)
  tas.analysis() 

  mseTotal = mseClass(tas.post, huss.post) 
  return pr.post, mseTotal 






def mse_PlottingAnalysis(path, location, scenario, calLength):

  """
  location: list object containing concluding with "_loc" which has been pre-defined
            with lat and lon values in accordance with cdo requirements for re-mapping 
  scenario: String with experiment specifics (4x, 1%, etc.)          

  ### Trouble with getting "site"_pr to execute, so caved and programmed as following 
  """

  prTotal, mseTotal = makeMSE(path, location, scenario, calLength)

  for site in location.keys(): 
    print "Now processing " + str(site)

    sitename = site.split("_")[0]

    # -- get location values for pr and MSE 
    mseTotal.mse_pixel(location[site])
    sitename_mse = mseTotal.pixel  
     
    lonSize = 2.5
    latSize = 2  

    # -- take larger grid box range in lieu of individual point  
    thisLat = float(location[site].split("lat=")[1])
    thisLon = float(location[site].split("lon=")[1].split("_lat")[0]) 
    areaRange = str(thisLon - float(lonSize))+","+str(thisLon + float(lonSize))+","+str(thisLat - float(latSize))+","+str(thisLat + float(latSize))
    
    sitename_pr = cdo.fldmean(input = cdo.sellonlatbox(areaRange, input = prTotal))  
    #sitename_pr = cdo.remapnn(location[site], input = pr.post)
    multiContourPlot(sitename, sitename_pr, sitename_mse, 15, ["Decadal Average ", "MSE ", "for ", titlecase(site), " Under ", str(scenario)], calLength)
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


graphics_path = os.path.join("/Users/xadx/Desktop/Dropbox", "MSE_Monsoon", "Area_Avg")


# -- GISS models
giss_4x = os.path.join(basePath, "abrupt4xCO2:E2-H_abrupt4xCO2_r1i1p3_day")
giss_rcp85 = os.path.join(basePath, "rcp85:E2-H_rcp85_r2i1p1_day")
giss_e2h_hist = os.path.join(basePath, "E2-H_historical_r6i1p1_day")
giss_e2r_hist  = os.path.join(basePath, "E2-R_historical_r6i1p1_day")

# -- Hadley models
had_base = "/Users/xadx/Desktop/Active_Research/Climate_Modeling/CMIP5_Data/HadGEM2"
had_1pct = os.path.join(had_base, "1pctCO2", "HadGEM2-ES_1pctCO2_r1i1p1_day")
had_4x = os.path.join(had_base, "abrupt4xCO2", "HadGEM2-ES.abrupt4xCO2.day.atmos.day.r1i1p1.v20130214")


hadrcp_dict = {} 
# -- multiple rcp 8.5 runs 
had_rcp85 = dict() 
for x in range(1, 5):
  print str(x)
  had_rcp85[x] = os.path.join(had_base, "rcp85", "HadGEM2-ES_rcp85_r" + str(x) + "i1p1")
  hadrcp_dict['had_rcp85'+str(x)]= [had_rcp85[x], "RCP 8.5 (HadGEM2) Run " + str(x)]

hadexp_dict = {'had1pct':[had_1pct, "1% CO2 (HadGEM2)"], 
      'had4x':[had_4x, "Abrupt 4x CO2 (HadGEM2)"]}

gissexp_dict = {'giss4x':[giss_4x, "Abrupt 4x CO2 (GISS)"], 
      'gissRCP85':[giss_rcp85, "RCP 8.5 (GISS)"],
      'e2hhist':[giss_e2h_hist, "Historical (GISS E2-H)"],
      'e2rhist':[giss_e2r_hist, "Historical (GISS E2-R)"]}


# -- coordinates of interest
mumbai_loc = "lon=73.01_lat=19.03"
nagpur_loc = "lon=79.0882_lat=21.1458"
kolkata_loc = "lon=88.3630400_lat=22.5626300" 

location_dict = {'mumbai':mumbai_loc, 
      'nagpur':nagpur_loc,
      'kolkata':kolkata_loc}




# -- at current moment am not distinguishing between grid box sizes of GISS v. Hadley models 


"""
# -- generate coordinate inputs for a spatial average over a geographic subset 
location_sellon_dict = {} 

for j in location_dict.keys(): 
  thisLat = float(location_dict[j].split("lat=")[1])
  thisLon = float(location_dict[j].split("lon=")[1].split("_lat")[0]) 
  thisLatRange = [thisLat - float(latSize), thisLat + float(latSize)]
  thisLonRange = [thisLon - float(lonSize), thisLon + float(lonSize)]
  location_sellon_dict[j] = str(thisLon - float(lonSize))+","+str(thisLon + float(lonSize))+","+str(thisLat - float(latSize))+","+str(thisLat + float(latSize))
 # coords = thisLonRange, thisLatRange



#kolkata_avg = cdo.fldavg(input = cdo.sellonlatbox(location_sellon_dict['kolkata'], input = "/var/folders/tl/3sdxh9rj0zd88dbwzj8bfjgc0000gn/T/cdoPyF7eRpd"))

"""



"""

To select the region with the longitudes from 120E to 90W and latitudes from 20N to 20S from all
input elds use:

cdo sellonlatbox,120,-90,20,-20 ifile ofile

"""


# -- location and experiment dictionaries 


hadcombined_dict = dict(hadexp_dict.items() + hadrcp_dict.items())

# -- ensure that all directories for figures exist, else create 
z = dict(hadcombined_dict.items() + gissexp_dict.items())
for i in z: 
  exp_path = os.path.join(graphics_path,z[i][1])
  if not os.path.exists(exp_path):
    os.makedirs(exp_path)



# -- each experiment-location pair for Hadley models
for exp in hadcombined_dict.keys():
  mse_PlottingAnalysis(hadcombined_dict[exp][0], location_dict, hadcombined_dict[exp][1], 360)


# -- each experiment-location pair for GISS models 
for exp in gissexp_dict.keys():
  mse_PlottingAnalysis(gissexp_dict[exp][0], location_dict, gissexp_dict[exp][1], 365)







# -- built as sample set to construct boxplot from 
prKolkata, mseKolkata = makeMSE(gissexp_dict['giss4x'][0], location_dict, gissexp_dict['giss4x'][1], 365)



df1 = pd.DataFrame(pctile, columns = ['Cume Rainfall'])
df1['DOY'] = range(364)
df1.plot(kind = 'scatter', x = 'Cume Rainfall', y = 'DOY')



# -- convert labels to months if either 360 or 365 day calendar : include if statement in there
months = [calendar.month_name[i] for i in range(1,13)]





# -- follow this procedure for editing labels to month values 
fig, ax = plt.subplots()
fig.canvas.draw() 
a=ax.get_xticks().tolist()
anew = [int(i*100+1) for i in a] 
adate = [datetime.strptime(str(x), '%j').strftime('%B %d') for x in anew]
ax.set_xticklabels(adate)
plt.show() 


# -- modify label markings to get month value 
datetime.strptime(x, '%j').strftime('%B %d')


range(1,13)
firstDay = [str(x)+"/"+str(1) for x in range(1,13)]
doyVals = [int(datetime.strptime(x, '%m/%d').strftime('%j')) for x in firstDay]
monthVals = [datetime.strptime(x, '%m/%d').strftime('%B') for x in firstDay] 



prKOLK = cdo.remapnn(kolkata_loc, input = prKolkata)


def makeBoxPlot(pr_series, Years, ylevel):
  # -- currently based on a single year, will have to develop it around 
  # -- taking decadal average
  # years: a list of (potentially) multiple year values 
  # returns a boxplot object that will be included inside an MSE plot
  # pr_series: CDO object of precipitation -- not making similar plots for any other series at the current moment   
  # ylevel: manual input to fix boxplot inside another plot 
  # currently setup as developing boxplot for single year -- will need to be modified in the future 

  pr_nc = Dataset(pr_series)
  precip = pr_nc.variables['pr'][:]
  #prK_time = cdo.showdate(input = pr_series)[0].split()
  pr_time = [datetime.strptime(x, '%Y-%m-%d') for x in cdo.showdate(input = pr_series)[0].split()]

  ct = 0 
  for y in Years:
    # -- index values corresponding to beginning and end of calendar year 
    first_val = [i.year == y for i in pr_time].index(True)
    last_val = [i for i, j in enumerate(pr_time) if j.year == y][-1] 

    # single year in vector format 
    thisYear = np.squeeze(precip[first_val:last_val])

    # -- statistics on annual total 
    total = sum(thisYear)
    pctile = np.cumsum(thisYear)/total 

    # -- creating boxplot from these parameters
    pct_val = [0.05, 0.25, 0.5, 0.75, 0.95]
    pct_val_dict = {}

    box_wider_big = np.array([-1,1,1,-1,-1])
    box_wider_sm = np.array([-1,1])

    for j in pct_val: 
      pct = [i > j for i in pctile]
      pct_val_dict[j] =  [l == True for l in pct].index(True)   

    # -- compute and plot specifics of the year's distribution   
    Mean=pct_val_dict[0.5]#mean
    IQR=[pct_val_dict[0.25],pct_val_dict[0.75]] #inter quantile range
    CL=[pct_val_dict[0.05],pct_val_dict[0.95]] #confidence limit
    A=np.random.random(50)
    D=plt.boxplot(A, vert = False) # a simple case with just one variable to boxplot
    D['medians'][0].set_xdata(pct_val_dict[0.5])
    # medians is 2x2
    D['medians'][0]._xy[:,1]=ylevel[ct]+box_wider_sm
   # D['medians']
   # boxes is 5rowx2
    D['boxes'][0]._xy[[0,1,4],0]=IQR[0]
    D['boxes'][0]._xy[[2,3],0]=IQR[1]
    D['boxes'][0]._xy[:,1]=ylevel[ct]+box_wider_big
    print D['boxes'][0]._xy 
    # whiskers is 2x2
    D['whiskers'][0].set_xdata(np.array([IQR[0], CL[0]]))
    D['whiskers'][0]._xy[:,1]=ylevel[ct]
    D['whiskers'][1].set_xdata(np.array([IQR[1], CL[1]]))
    D['whiskers'][1]._xy[:,1]=ylevel[ct]
    print D['whiskers'][0]
    # caps is 2x2
    D['caps'][0].set_xdata(np.array([CL[0], CL[0]]))
    D['caps'][0]._xy[:,1]=ylevel[ct]+box_wider_sm
    D['caps'][1].set_xdata(np.array([CL[1], CL[1]]))
    D['caps'][1]._xy[:,1]=ylevel[ct]+box_wider_sm




    # -- add text of years to plot 
    plt.text(pct_val_dict[0.5]-5,ylevel[ct]-2, str(y)+"-"+str(y+10), fontsize = 16)

    ct += 1 # counter for identifying where the box plot should be vertically placed 
    _=plt.xlim(np.array(CL)+[-0.1*np.ptp(CL), 0.1*np.ptp(CL)]) #reset the limit
  
  # -- modify x-axis labels 
  middleOfMonth = [str(x)+"/"+str(15) for x in range(1,13)]
  doyVals = [int(datetime.strptime(x, '%m/%d').strftime('%j')) for x in middleOfMonth]
  monthVals = [datetime.strptime(x, '%m/%d').strftime('%B %d') for x in middleOfMonth] 

  plt.xticks(doyVals, monthVals)
  _=plt.xlim(np.array(CL)+[-0.1*np.ptp(CL), 0.1*np.ptp(CL)]) #reset the limit
  _=plt.ylim(min(ylevel)-2.5, max(ylevel)+2.25)
  return D 


kol_bp = makeBoxPlot(prKOLK, [2900, 2960], [330,320])
plt.show() 






def decLinePlotwBox(in_data, titleText, subset, calLength, pr_data, pr_year, verticalShift):
  """

  'indata' is a single series in a CDO format, - either from a single pixel
            or spatial average of multiple pixels  
  either from a zonal average or individual location 
  'titleText' is a character string that labels plot title 
  calLength = either 360 or 365 depending on model used 

  pr_data is the precip series that corresponds to the series from which the MSE data is constructed
  pr_year is the individual year from which the boxplot gets constructed
  """ 
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

  # -- need range to get full interval, not just end points
  decmean = [z[:,range(indexes[:,i][0],indexes[:,i][1])].mean(axis=1) for i in range(indexes.shape[1])] 
  plt.subplots_adjust(hspace=0.07)

  pnum = 1 
  f = 1 
  x_pr_restrict = [150,275]  # portion kept in final plots 
  x_restrict = [0,calLength]

  # -- take average across individual decades, and label accordingly 
  for c in range(indexes.shape[1]):
    print "Now processing " + str(c) + " of the decadal mean."

    znew[:,c] = decmean[c]

    # -- labeling of correct decades, intervals are decades  
    year1 = int(years[0]) + indexes[0,c] -1
    year2 = int(years[0]) + indexes[1,c] -1  

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

  makeBoxPlot(pr_data, pr_year, verticalShift)   
  #plt.tight_layout()  
  fn = os.path.join(graphics_path, titleText[5], titleText[1] + "_" + titleText[3] + "_" + titleText[5] + "_Decadal_Plot.eps")
  plt.savefig(fn, figsize=(8, 6), dpi = 1200, transparent = True, facecolor='w', edgecolor='k') 
  plt.close()  
  #plt.show()


### TESTING CODE ABOVE 


# first get a single location, experiment
rcp85PR, rcp85MSE = makeMSE(giss_rcp85, mumbai_loc, "RCP 8.5 (GISS)", 365)
# then get the single remappedNN series 
rcp85PR_mumbai = cdo.remapnn(mumbai_loc, input = rcp85PR)
rcp85MSE_mumbai = cdo.remapnn(mumbai_loc, input = rcp85MSE)

rcp85MSE.mse_pixel(mumbai_loc)

# find 


makeBoxPlot(rcp85PR_mumbai, 2015, 330)
decLinePlotwBox(rcp85MSE.pixel, ["Decadal Average ",  "MSE ", "for ", "Mumbai", " Under ", "RCP 8.5"], False, 365, rcp85PR_mumbai, 2015, 330)










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
