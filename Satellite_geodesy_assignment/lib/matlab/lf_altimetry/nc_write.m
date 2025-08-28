%matlab script to write a two dimensional grid z to a file
% adapted from fortran subroutine NC_write
%Roelof Rietbroek 02-10-2007
%usage
% nc_write(file,x,y,z,varargin)
% file is the filename
% x longitude
% y latitude
% z values z(nlon,nlat)
% varargin additional parameters
%

%adapted 13-09-2010
%use new netcdf interface
function nc_write(file,x,y,z,varargin)

%defaults



history='unknown';

X.Name='x';
X.Nctype='NC_DOUBLE';
X.Dimension={'x'};
X.Attribute(1).Name='units';
X.Attribute(1).Value='degrees east';
X.Attribute(2).Name='long_name';
X.Attribute(2).Value='Longitude';
X.Attribute(3).Name='actual_range';
X.Attribute(3).Value=[min(x) max(x)];

Y.Name='y';
Y.Nctype='NC_DOUBLE';
Y.Dimension={'y'};
Y.Attribute(1).Name='units';
Y.Attribute(1).Value='degrees north';
Y.Attribute(2).Name='long_name';
Y.Attribute(2).Value='Latitude';
Y.Attribute(3).Name='actual_range';
Y.Attribute(3).Value=[min(y) max(y)];


Z.Name='z';
Z.Nctype='NC_DOUBLE';
Z.Dimension={'y','x'};
Z.Attribute(1).Name='units';
Z.Attribute(1).Value='unknown';
Z.Attribute(2).Name='long_name';
Z.Attribute(2).Value='unknown';
Z.Attribute(3).Name='actual_range';
Z.Attribute(3).Value=[min(min(z)) max(max(z))];


ctime=[];
%process variable arguments

nargs=nargin-4;
if(nargs ~= 0)
  for i=1:nargs
    switch varargin{i}
     case 'history',
      history=varargin{i+1};
     case 'zunits',
      Z.Attribute(1).Value=varargin{i+1};
     case 'xunits',
      X.Attribute(1).Value=varargin{i+1};
     case 'yunits',
      Y.Attribute(1).Value=varargin{i+1};
     case 'xname',
      X.Attribute(2).Value=varargin{i+1};
     case 'yname',
       Y.Attribute(2).Value=varargin{i+1};
     case 'zname',
       Z.Attribute(2).Value=varargin{i+1};

     case 'ctime',
      ctime=varargin{i+1};
      
    end
  end
end

%create new file
ncid = netcdf.create(file,'CLOBBER');


%set global attributes
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history',history);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','COARDS');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'node_offset',0);
if(not(isempty(ctime)))
  netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Time_year',ctime);
end

%setup dimensions (x,y)
xdid=netcdf.defDim(ncid,'x', length(x));
ydid=netcdf.defDim(ncid,'y', length(y));


%define variables
xid = netcdf.defVar(ncid,X.Name,X.Nctype,xdid);
for i=1:3
  netcdf.putAtt(ncid,xid,X.Attribute(i).Name,X.Attribute(i).Value);
end


yid = netcdf.defVar(ncid,Y.Name,Y.Nctype,ydid);
for i=1:3
  netcdf.putAtt(ncid,yid,Y.Attribute(i).Name,Y.Attribute(i).Value);
end

zid = netcdf.defVar(ncid,Z.Name,Z.Nctype,[xdid ydid]);
for i=1:3
  netcdf.putAtt(ncid,zid,Z.Attribute(i).Name,Z.Attribute(i).Value);
end

netcdf.endDef(ncid);

%add data to netcdf grid
netcdf.putVar(ncid,xid,x);
netcdf.putVar(ncid,yid,y);
netcdf.putVar(ncid,zid,z);

netcdf.close(ncid);



