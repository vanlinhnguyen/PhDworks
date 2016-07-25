function ncdump(theNetCDFFile)

nc = netcdf(theNetCDFFile);%open for reading

%disp all informations (dimensions, variables, attributes, ...) of the netcdf file
nc.getInfo

%close nc
close(nc);