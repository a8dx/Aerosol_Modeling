# -- cdoCl.py 



import glob 
from cdo import * 

cdo = Cdo()



class cdoClass:

  

  def __init__(self, var, IndiaSubset):
    self.var = var
    self.series = var+"_day"
    self.files = []
    self.subset = IndiaSubset

# -- searches path and identifies relevant files for specified var name
  def get_files(self): 
    self.files =  sorted(glob.glob(self.series + '_*.nc'))
    print "Files in this segment are " + str(self.files)
    print "Files sorted? " + str(sorted(self.files) == self.files)

# -- concatenates files across all years 
  def cat_files(self):
    print "\n Beginning concatenation. \n " 
    self.catfiles = cdo.cat(input = self.files)   

  def india_subset(self):
    print "Subset to India bounding box alone? " + str(self.subset)
    if (self.subset == True):
      bbox = "53.75,93.75,1,41"
      print "\n Restricting coordinates to " + str(bbox) + "\n"     
      self.post = cdo.sellonlatbox(str(bbox), input = self.catfiles)
    else: 
      print "\n Maintaining original coordinates of " + str(self.files) + "\n"
      self.post = self.catfiles

# -- completes the analysis by invoking all indicated subroutines 
  def analysis(self):
    print "\n Conducting full analysis on " + str(self.var) + " series. \n " 
    self.get_files()
    self.cat_files() 
    self.india_subset() 
    print "\n ***** Analysis has been successfully completed ***** \n"

