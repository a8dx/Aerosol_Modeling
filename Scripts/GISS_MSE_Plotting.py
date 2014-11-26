# Plotting functions for use by GISS_MSE_Analysis.py 


from cdo import * 
cdo = Cdo()


def contourPlot(in_data, num_levels, titleText):
  """

  'indata' is a single series in a CDO format, 
  either from a zonal average or individual location 
  'num_levels' is an integer for the number of levels in the contour plot
  'titleText' is a character string that labels plot title 
  """ 
  
  # -- convert years into integers
  indata = Dataset(in_data)
  yearStr = cdo.showyear(input = in_data)[0].split() 
  yearRange = [int(i) for i in yearStr]
  time_vals = (1461.0/1460.0)*indata.variables['time'][:]
  time_dateFormat = [relativedelta(years = int(yearRange[0]) - 1) + datetime.date.fromordinal(int(round(i))) for i in time_vals]   # convert into date-time format 
 # labels = [i.strftime("%B")[:3] + "-" + i.strftime("%d") for i in time_dateFormat]
  labels = [i.isoformat() for i in time_dateFormat] 

  # -- contour data series
  z = indata.variables[cdo.showname(input = in_data)[0]][:]
  zf = z.squeeze().reshape(len(time_vals)/len(yearRange),len(yearRange), order = 'F')

  # -- done to ensure the calculations are correct
  y = time_vals.reshape(len(yearRange), len(time_vals)/len(yearRange))
  dayRange = y[0,:]
  X,Y = np.meshgrid(yearRange, dayRange)

  day_dateFormat = [relativedelta(years = int(yearRange[0]) - 1) + datetime.date.fromordinal(int(round(i))) for i in dayRange]
#  dayLabels = [i.strftime("%B")[:3] + "-" + i.strftime("%d") for i in day_dateFormat]
  dayLabels = [j.isoformat() for j in day_dateFormat] 


  #pl.axes([0.05, 0.05, 0.90, 0.90])

  fig, ax = plt.subplots(1,1)
  plt.xlabel("Year")
  plt.ylabel("Day")
  plt.title(titleText)

  #levels = np.linspace(310,370,7)

  cs = ax.contourf(X,Y,zf,
    num_levels, 
    #levels = levels, 
    cmap = plt.cm.hot_r)
  fig.colorbar(cs, ax=ax, format = "%.3f")
  plt.show() 
 






def multiContourPlot(location, var1, var2, title1, title2, scen):
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
  title1/2: Titles correspond to variables being analyzed, original sequence is 
            precip for slot 1 and MSE for slot 2 
  scen: string featuring experiment/scenario (4xCO2, 1%CO2, etc.)           

  """ 

  bnds = [5,10,20] 
  [contourPlot(var1, x, location.capitalize() + str(title1) + str(scen)) for x in bnds]
  [contourPlot(var2, x, location.capitalize() + str(title2) + str(scen)) for x in bnds]

#  [contourPlot(eval(str(location)+"_"+str(var1)), x, location.capitalize() + str(title1) + str(scen)) for x in bnds]
#  [contourPlot(eval(str(location)+"_"+str(var2)), x, location.capitalize() + str(title2) + str(scen)) for x in bnds]






def decLinePlot(in_data, titleText):
  """

  'indata' is a single series in a CDO format, - either from a single pixel
            or spatial average of multiple pixels  
  either from a zonal average or individual location 
  'titleText' is a character string that labels plot title 
  """ 
  
    # -- convert years into integers
  indata = Dataset(in_data)
  yearStr = cdo.showyear(input = in_data)[0].split() 
  yearRange = [int(i) for i in yearStr]
  time_vals = (1461.0/1460.0)*indata.variables['time'][:]
  time_dateFormat = [relativedelta(years = int(yearRange[0]) - 1) + datetime.date.fromordinal(int(round(i))) for i in time_vals]   # convert into date-time format 
  # labels = [i.strftime("%B")[:3] + "-" + i.strftime("%d") for i in time_dateFormat]
  labels = [i.isoformat() for i in time_dateFormat] 

  # -- contour data series
  z = indata.variables[cdo.showname(input = in_data)[0]][:]
  zf = z.squeeze().reshape(len(time_vals)/len(yearRange),len(yearRange), order = 'F')

  # -- done to ensure the calculations are correct
  y = time_vals.reshape(len(yearRange), len(time_vals)/len(yearRange))
  dayRange = y[0,:]
  X,Y = np.meshgrid(yearRange, dayRange)

  # -- total number of decadal periods 
  it = math.floor(zf.shape[1]/10)
  indexes = np.array([range(int(it)), range(1,1+int(it))])*10



  # -- initialize output array 
  znew = np.zeros((zf.shape[0], int(it)))


  m = ["o", "s", "o", "s", "o"]

  plt.subplots(1,1) 

  decmean = [zf[:,indexes[:,i]].mean(axis=1) for i in range(indexes.shape[1])] 

  # -- take average across individual decades, and label accordingly 
  for c in range(indexes.shape[1]):
    znew[:,c] = decmean[c]
    #ax = fig.add_subplot(111) 
    plt.plot(dayRange, znew[:,c], m[c])


  # -- need to limit the y-extent of precip.  
  plt.show() 

  # -- include something about labels here 






def eof_standard_plot(obj, value, var, titleText):
  """ Produces some standard plot output and therefore easier to change 
      these settings than manual work. 
      Figure out some way to default to total range for index (i.e., not just 1)
      obj : post-processed netCDF4 object (i.e., not raw original object)
      var : prescribed abbreviation of climate variable of interest 
      value : when working with array, take a single value (i.e., 1st EOF)

      NEED TO INCLUDE FILENAME, CURRENTLY COMMENTED OUT
  """
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  cliva = obj.variables[var][value,:,:].squeeze()   # climate variable
  cliva_units = obj.variables[var].units

  

  lon_0 = lons.mean()
  lat_0 = lats.mean()

  m = Basemap(width=5000000,height=3500000,
            resolution='l',projection='stere',\
            lat_ts =40, lat_0=lat_0,lon_0=lon_0) 

  # Because our lon and lat variables are 1D, 
  # use meshgrid to create 2D arrays 
  # Not necessary if coordinates are already in 2D arrays.
  lon, lat = np.meshgrid(lons, lats)
  xi, yi = m(lon, lat)

  # Plot Data
  cs = m.pcolormesh(lon,lat,cliva,shading='flat',cmap=plt.cm.jet,latlon=True)

  # Add Grid Lines
  m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
  m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

  # Add Coastlines, States, and Country Boundaries
  m.drawcoastlines()
  m.drawstates()
  m.drawcountries()

  # Add Colorbar
  cbar = m.colorbar(cs, location='bottom', pad="10%")
  xu = str(cliva_units) + " yr-1"
  cbar.set_label(xu)

  # Add Title
  plt.title(titleText)
  #plt.show()
  #plt.draw()   
  #fn = filename + ".eps"     # hard-coding here can be improved
  #plt.savefig(fn)
  #plt.show()
  #fig.show() 
  #def ensemble_avg(modphys)
  show() 
  obj.close()  # closes the file




def kav7Plot(obj, value, var, titleText):
  """ Produces some standard plot output and therefore easier to change 
      these settings than manual work. 
      Figure out some way to default to total range for index (i.e., not just 1)
      obj : post-processed netCDF4 object (i.e., not raw original object)
      var : prescribed abbreviation of climate variable of interest 
      value : when working with array, take a single value (i.e., 1st EOF)

      NEED TO INCLUDE FILENAME, CURRENTLY COMMENTED OUT
  """
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  cliva = obj.variables[var][value,:,:].squeeze()   # climate variable
  cliva_units = obj.variables[var].units
  lon_0 = lons.mean()

  m = Basemap(projection = 'kav7', lon_0 = lon_0, resolution = None)
  m.drawmapboundary(fill_color='0.3')

  # Because our lon and lat variables are 1D, 
  # use meshgrid to create 2D arrays 
  # Not necessary if coordinates are already in 2D arrays.
  lon, lat = np.meshgrid(lons, lats)
  xi, yi = m(lon, lat)

  # Plot Data
  cs = m.pcolormesh(lon,lat,cliva,shading='flat',cmap=plt.cm.jet,latlon=True)

  # Add Grid Lines
  m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
  m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

  # Add Colorbar
  cbar = m.colorbar(cs, location='bottom', pad="10%")
  cbar.set_label(cliva_units)

  # Add Title
  plt.title(titleText)
  show() 
  obj.close()  # closes the file






def simplePlot(obj, var, titleLabel):
  """
    Highly convertible definition to produce plots on the fly from netCDF4 objects 
    and user-specified climate variables. 

    Assumes that yvar is in a (N,) array format (i.e., 1D object) 
    obj:  post-'Dataset' object (until figuring out how to Dataset inside function)
    time: variable plotted on x-axis
    yvar: variable plotted on y-axis
    lab:  plot label text  
  """

  yvar = var
  #lab = "Insert text here"
  time = obj.variables['time'][:]
  dt_time = [dt.date(1850,1,1) + dt.timedelta(days=t) for t in time]
  #time.shape = time.shape[0]    # keep first dimension only
  yt = obj.variables[yvar][:]
  yt.shape = yt.shape[0]    # keep first dimension only
  yt_units = obj.variables[yvar].units
  #ticks = num2date(t1, "days since 1850-01-01 00:00:00")
  #fig = figure()
  title(titleLabel)
  plot(dt_time,yt, label = title)
  xlabel("Time")
  ylabel(yt_units)
  xticks(rotation = 55)
  #autofmt_xdate()
  #show()
  #obj.close()


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

#  count = 0   # not certain what role of this is   
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
 #   count += 1    # identifier for ensemble total 

    if (plotMembersClima == True):
      ensemblePlot()

    print str(len(ensembleTotal)) + " length of ensemble total"  
  print ensembleTotal  
  return ensembleTotal





def monthlyPlotDeprecated(object, var, title):
  
  # time goes 0:11

  obj = Dataset(object)
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  cliva_units = obj.variables[var].units

  lon_0 = lons.mean()
  lat_0 = lats.mean()
  

  fig, axes = plt.subplots(nrows = 4, ncols = 3)

  count = 0    # to identify slices 
  for ax in axes.flat:
    
    cliva = obj.variables[var][count,:,:]   # climate variable

    #m_ax = Basemap(ax = ax, width=5000000,height=3500000,resolution='l',projection='stere', lat_ts=40,lat_0=lat_0,lon_0=lon_0)
    m_ax = Basemap(ax = ax, width=5000000,height=3500000, projection = 'stere', lat_ts = 40, lat_0 = lat_0, lon_0 = lon_0)

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
 


  # grid = AxesGrid(fig, 132, # similar to subplot(132)
  #                   nrows_ncols = (2, 2),
  #                   axes_pad = 0.0,
  #                   share_all=True,
  #                   label_mode = "L",
  #                   cbar_location = "top",
  #                   cbar_mode="single",
  #                   )

  # for i in range(12):
  #   im = grid[i].imshow(Z, extent=extent, interpolation="nearest")
  #   #plt.colorbar(im, cax = grid.cbar_axes[0])
  # grid.cbar_axes[0].colorbar(im)

  # for cax in grid.cbar_axes:
  #   cax.toggle_label(False)
        
  #   # This affects all axes as share_all = True.
  # grid.axes_llc.set_xticks([-2, 0, 2])
  # grid.axes_llc.set_yticks([-2, 0, 2])

  plt.suptitle(title)
  plt.tight_layout()
  cax = plt.axes([0.1, 0.05, 0.8, 0.025])
  plt.colorbar(cax=cax, orientation = 'horizontal')  
  plt.subplots_adjust(left = 0.1, right = 0.9, top = 0.1, bottom = 0.9)
  plt.show()   

def monthlyPlot(object, var, title):
  
  obj = Dataset(object)
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  times = obj.variables['time'][:]
  cliva_units = obj.variables[var].units
  lon_0 = lons.mean()
  lat_0 = lats.mean()
  clevs = np.arange(-30,30.1,2.)
  lon, lat = np.meshgrid(lons, lats)  
  m = Basemap(width=4000000,height=3500000,resolution='l',projection = 'stere', lat_ts = 5, lat_0 = lat_0, lon_0 = lon_0)  #  
  xi, yi = m(lon, lat)
  fig=plt.figure()

  for nt,time in enumerate(times):
    ax = fig.add_subplot(4,3,1+nt)
    cs = m.pcolor(xi,yi,np.squeeze(obj.variables[var][nt,:,:]))
    #cs = m.contourf(xi,yi,np.squeeze(obj.variables[var][nt,:,:]),clevs,cmap=plt.cm.jet,extend='both')
    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

    plt.title(calendar.month_abbr[nt+1],fontsize=9)

  plt.figtext(0.5,0.95,title,horizontalalignment='center',fontsize=14)
  #plt.tight_layout()
  cax = plt.axes([0.1, 0.05, 0.8, 0.025])
  cbar = plt.colorbar(cax=cax, orientation='horizontal')
  cbar.set_label(cliva_units)
  cbar.ax.set_xlabel(cliva_units)
  plt.show()   



def rectMonthlyPlot(object, var, title):
  
  obj = Dataset(object)
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  times = obj.variables['time'][:]
  cliva_units = obj.variables[var].units
  lon_0 = lons.mean()
  lat_0 = lats.mean()
  clevs = np.arange(-30,30.1,2.)
  lon, lat = np.meshgrid(lons, lats)  
  m = Basemap(width=4000000,height=3500000,resolution='l',projection = 'stere', lat_ts = 5, lat_0 = lat_0, lon_0 = lon_0)  #  
  xi, yi = m(lon, lat)
  fig=plt.figure()

  for nt,time in enumerate(times):
    ax = fig.add_subplot(4,3,1+nt)
    cs = m.pcolor(xi,yi,np.squeeze(obj.variables[var][nt,:,:]))
    #cs = m.contourf(xi,yi,np.squeeze(obj.variables[var][nt,:,:]),clevs,cmap=plt.cm.jet,extend='both')
    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(np.arange(-80., 81., 10.), labels=[1,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1], fontsize=10)

    plt.title(calendar.month_abbr[nt+1],fontsize=9)

  plt.figtext(0.5,0.95,title,horizontalalignment='center',fontsize=14)
  #plt.tight_layout()
  cax = plt.axes([0.1, 0.05, 0.8, 0.025])
  cbar = plt.colorbar(cax=cax, orientation='horizontal')
  cbar.set_label(cliva_units)
  cbar.ax.set_xlabel('UNITS HERE')
  plt.show()   



def spacePlot(obj, slice, var):
  """
  Will assume post-Dataset for now, until can get working
  Currently not saving plot to file - uncomment to do so.   
  """
  lons = obj.variables['lon'][:]
  lats = obj.variables['lat'][:]
  cliva = obj.variables[var][0,:,:]   # climate variable
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
  plt.show() 
  # Add Title
  #plt.title(titleText)
  #plt.show()
  #plt.draw()   
  #fn = "EOF1" + str(goodFiles[0][0:-3]) + ".eps"     # hard-coding here can be improved
  #plt.savefig(fn)
  #plt.show()
  #fig.show() 
  #def ensemble_avg(modphys)
 # obj.close()  # closes the file  



def ensemblePlot():
  """
  ensemblePlot: Appends monthly climatologies for each ensemble member 
  as looping through a vector with concatenated member files. 
  Currently requires no input since it feeds off output objects from 
  standard `genEnsembleAvg' function.
  """ 
  ensNumber = 1    # extract ensemble member number from the file path eventually 
  for member in ensembleTotal:
    y = cdo.ymonavg(input = "-fldavg " + member) #, options = '-f nc')
    df = Dataset(y, mode = 'r')  
    dfTime = df.variables['time'][:]    #climatology time 
    dfVar = df.variables[varSeries][:]
    dfVar.shape = (12)   # reshape to 1D 
    seriesLab = str(modelPhys) + "Member " + str(ensNumber)
    plot(dfT, dfVar, label = seriesLab)
    ensNumber += 1 
    legend(loc = 'upper left')
    show() 








