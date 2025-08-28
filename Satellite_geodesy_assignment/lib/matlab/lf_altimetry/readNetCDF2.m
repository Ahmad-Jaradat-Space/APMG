function data=readNetCDF2(filename)
% Reads netCDF Version 3 and 4 files in Matlab
%
% INPUT: filename ... filename including path
%
% OUTPUT: data ... struct with netCDF data
%
% Note: All scale factors, add_offsets etc. will be applied so the
% resulting data should have the units decribed in the "units" cell array
% which is included inside the data struct.
%
% Uebbing 06/2016


if isOctave  %octave
    pkg load netcdf
    ncid = netcdf_open(filename, 'NC_NOWRITE'); % oeffnen der Datei
    [ndims,nvars,ngatts,unlimdimid] = netcdf_inq(ncid);
    for j=0:nvars-1
        [varname,xtype,dimids,natts] = netcdf_inqVar(ncid,j);
        data.(varname) = double(netcdf_getVar(ncid, j)); % einlesen der daten
        
        scaleFac = 1;
        addOff=0;
        fillVal=NaN;
        for k=0:natts-1
            attname = netcdf_inqAttName(ncid,j,k);
            switch attname
                
                case 'scale_factor'
                    scaleFac = netcdf_getAtt(ncid,j,'scale_factor');  % Scaling factor
                case 'add_offset'
                    addOff = netcdf_getAtt(ncid,j,'add_offset'); % Additionskonstante
                case '_FillValue'
                    fillVal = netcdf_getAtt(ncid,j,'_FillValue'); % Fillvalue (NaN)
                case 'units'
                    units{j+1} = netcdf_getAtt(ncid,j,'units'); % Einheit
                otherwise
                    %nix
                    
            end
            
        end
        
        % Transformieren der Daten in richtige Einheit
        if (isnan(fillVal)==0)
            data.(varname)(data.(varname)==fillVal)=NaN;
        end
        data.(varname) = data.(varname)*double(scaleFac) + double(addOff);
        
        %     data.(varname)=adat;
        
    end
    
    try
        data.units = units;
    catch
        
    end
    try
        % read global attributes
        for j=0:ngatts-1
            attname = netcdf_inqAttName(ncid,netcdf_getConstant('NC_GLOBAL'),j);
            [xtype, attlen] = netcdf_inqAtt(ncid,netcdf_getConstant('NC_GLOBAL'),attname);
            globAtt.(attname)=netcdf_getAtt(ncid,netcdf_getConstant('NC_GLOBAL'),attname);
        end
        data.globAtt = globAtt;
    catch
    end
    netcdf_close(ncid);
else % matlab
    vv = version;
    vk = strfind(vv,'R');
    vv = str2double(vv(vk+1:vk+4));
    if vv<=2014
        ncid = netcdf.open(filename, 0); % oeffnen der Datei
    else
        ncid = netcdf.open(filename, 'NC_NOWRITE'); % oeffnen der Datei
    end
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for j=0:nvars-1
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,j);
        data.(varname) = double(netcdf.getVar(ncid, j)); % einlesen der daten
        
        scaleFac = 1;
        addOff=0;
        fillVal=NaN;
        for k=0:natts-1
            attname = netcdf.inqAttName(ncid,j,k);
            switch attname
                
                case 'scale_factor'
                    scaleFac = netcdf.getAtt(ncid,j,'scale_factor');  % Scaling factor
                case 'add_offset'
                    addOff = netcdf.getAtt(ncid,j,'add_offset'); % Additionskonstante
                case '_FillValue'
                    fillVal = netcdf.getAtt(ncid,j,'_FillValue'); % Fillvalue (NaN)
                case 'units'
                    units{j+1} = netcdf.getAtt(ncid,j,'units'); % Einheit
                otherwise
                    %nix
                    
            end
            
        end
        
        % Transformieren der Daten in richtige Einheit
        if (isnan(fillVal)==0)
            data.(varname)(data.(varname)==fillVal)=NaN;
        end
        data.(varname) = data.(varname)*double(scaleFac) + double(addOff);
        
        %     data.(varname)=adat;
        
    end
    
    try
        data.units = units;
    catch
        
    end
    try
        % read global attributes
        for j=0:ngatts-1
            attname = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),j);
            [xtype, attlen] = netcdf.inqAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
            globAtt.(attname)=netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'),attname);
        end
        data.globAtt = globAtt;
    catch
    end
    
    netcdf.close(ncid);
end

end



%%
%% Return: true if the environment is Octave.
%%
function retval = isOctave
persistent cacheval;  % speeds up repeated calls

if isempty (cacheval)
    cacheval = (exist ('OCTAVE_VERSION', 'builtin') > 0);
end

retval = cacheval;
end
