#!/bin/bash

ct=0
for f in *.*
do
	echo "Processing $f ..."

done 

for d in *E2-H*313*/
do
	echo "Identifying only historical Misc directories, such as $d"
	cd $d 
	for j in pr_*.nc 
	do
		echo "File processing is $j"
		#/usr/local/bin/cdo info "$j" 
		#/usr/local/bin/cdo cat FILES  output.nc   # must include file extension for this to work  
		#/usr/local/bin/cdo sellonlatbox,58.75,96.25,5,41 output.nc output_subset.nc 


		#
		ct=$(($ct+1))
		echo "Count is $ct"
	done 	

done 


