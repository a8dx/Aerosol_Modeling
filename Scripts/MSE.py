
# -- MSE.py 


import glob 
from cdo import * 
import numpy as np
import matplotlib.pyplot as plt

cdo = Cdo()

class mseClass: 

  def __init__(self, tasData, hussData):
    #self.var = self
    self.tas = tasData
    self.huss = hussData 

  def calcMSE(self):
    print "\n Calculating MSE for the given dataset \n"
    cpT = cdo.mulc(1004, input = self.tas)
    Lvq = cdo.mulc(2.5e6, input = self.huss)
    h = cdo.add(input = " ".join([cpT,Lvq]))  # MSE
    self.mse = cdo.divc(1004, input = h)

# -- generates a bar plot for a given location-year for the specified series 
  def mse_barplot(self, location, year, titleheading):
    print "\n Generating single-year barplot \n"    
    self.calcMSE() 
    self.year = year 
    self.location = location 
    self.title = titleheading 
    self.local = cdo.remapnn(str(self.location), input = self.mse)
    dS = Dataset(cdo.selyear(self.year, input = self.local)) 
    time_vals = dS.variables['time'][:]
    year0 = int(cdo.showyear(input = self.mse)[0][0:4])
    time_dateFormat = [relativedelta(years = year0 - 1) + datetime.date.fromordinal(int(round(i))) for i in time_vals]   # convert into date-time format 
    labels = [i.strftime("%B")[:3] + "-" + i.strftime("%d") for i in time_dateFormat]
    y_var = dS.variables[cdo.showname(input = self.local)[0]][:]
    y_var.shape = y_var.shape[0]
    bin = time_dateFormat
    width = 0.2
    ax = plt.subplot(111)
    ax.bar(bin, y_var, width, color = 'r')
    ax.set_ylim(bottom = 300, top = 375)
    title(str(self.title))
    xticks(rotation = 55)
    plt.show()

  def mse_scatterplot(self, location, years, titleheading):
    
    # years 
    self.calcMSE()
    self.location = location 
    self.years = years
    self.title = titleheading 

    # -- generate years to do a time subset with CDO 
    rg = "" 
    yrSpan = ()
    for yr in range(self.years[0],self.years[1]+1):  # -- inclusive 
      rg = rg+","+str(yr)
      yrSpan = yr
    rg = rg[1:]  # -- remove leading comma  
    
    locAllYears = cdo.remapnn(str(self.location), input = cdo.selyear(rg, input = self.mse))
    year0 = int(cdo.showyear(input = self.mse)[0][0:4])

    ax = plt.subplot(111)

    # -- create loop and for each year produce a line plot 
    for yrLoop in range(self.years[0],self.years[1]+1):
      singleYear = Dataset(cdo.selyear(yrLoop, input = locAllYears))
      self.time_vals = (1461.0/1460.0)*singleYear.variables['time'][:]
    #  self.time_vals = singleYear.variables['time'][:]
     
      time_dateFormat = [relativedelta(years = year0 - 1) + datetime.date.fromordinal(int(round(i))) for i in self.time_vals]   # convert into date-time format 
      print str(time_dateFormat) + "These are time_dateFormat objects"
      labels = [i.strftime("%B")[:3] + "-" + i.strftime("%d") for i in time_dateFormat]
      print str(labels) + "These are labels"
      y_var = singleYear.variables[cdo.showname(input = locAllYears)[0]][:]
      y_var.shape = y_var.shape[0]
      plt.plot(time_dateFormat, y_var, label = yrLoop)

    title(str(self.title))
    xticks(rotation = 55)
    legend = plt.legend(loc = 'lower center', shadow = False)
    show()   


  def mse_pixel(self, location):
    self.calcMSE() 
    self.location = location 
    self.pixel = cdo.remapnn(str(self.location), input = self.mse)


  def mse_time_subset(self, location, years):
    self.location = location    
    self.years = years 
    self.mse_pixel(self.location)
    # -- generate years to do a time subset with CDO 
    rg = "" 
    yrSpan = ()
    for yr in range(self.years[0],self.years[1]+1):  # -- inclusive 
      rg = rg+","+str(yr)
      yrSpan = yr
    rg = rg[1:]  # -- remove leading comma  
    self.time_subset = cdo.selyear(rg, input = self.pixel)

    
