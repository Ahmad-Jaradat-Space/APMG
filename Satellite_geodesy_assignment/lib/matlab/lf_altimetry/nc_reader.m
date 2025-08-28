function RES = nc_reader(filename)

ncid = netcdf.open(filename, 'NC_NOWRITE');

finfo = ncinfo(filename);
RES=struct();

%% read all attributes
for i=1:length(finfo.Attributes)
    eval(['RES.',finfo.Attributes(i).Name,'=finfo.Attributes(i).Value;']);
end

%% read all variables
for i=1:length(finfo.Variables)
    
    A.data = ncread(filename,finfo.Variables(i).Name);
    A.attributes=finfo.Variables(i).Attributes;
    
    eval(['RES.' finfo.Variables(i).Name '=A.data;']);
    
    
end


netcdf.close(ncid);
