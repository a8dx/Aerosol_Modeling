cd('/Users/xadx/Desktop/Active_Research/Climate_Modeling/GISS_Downloads/Precip/E2-R_historicalMisc_r5i1p313')


d = 'pr_Amon_GISS-E2-R_historicalMisc_r5i1p313_185001-187512.nc'  


nc = netcdf.open(d, 'NC_NOWRITE');
[dimname, dimlen] = netcdf.inqDim(ncid,0); 