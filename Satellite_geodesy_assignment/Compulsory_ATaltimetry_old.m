    %
% Test the GPOD TALES results 
% usedin SCOOP and all reprt (maintain same 4 figures in output file jpg and verify statistics)
% by: L. Fenoglio, Buchhaupt
% 09.08.2019 adaped to test TUDaBo results for CS2 tst track TALES
% 16.09.2020 Encapsulated plotting and statistics calculation in their own
%            functions, made the program use configuration files and if
%            multiple configs are provided it will loop over them and
%            execute each one (M. G�rtner)
% 29.09.2020 Added functionality to limit the used Latitudes in the test
%            Added code to ensure, that the same points within a given
%            tolerance (default 1e-8) are used. This currently checks
%            samness via Latitude (M. G�rtner)

dbstop if error
close all
clearvars;
clc
addpath(genpath('lib/matlab/'));

% Read all provided configs in the Config folder
%[configFileList] = getConfigsInFolder([pwd,'/Config/AlongTrack/']);
[configFileList] = getConfigsInFolder([pwd,'/Config/AlongTrack/']);

% Loop over the configs
for ii=1:length(configFileList)

config=readConfig(configFileList{ii});

if ~isfield(config,'nameFile1')
    error('Config must contain nameFile1 to give a correct name to the results')
end
if ~isfield(config,'nameFile2')
    error('Config must contain nameFile2 to give a correct name to the results')
end
if ~isfield(config,'testfile1')
    error('Config must contain testfile1 to load the first file to test')
end
if ~isfield(config,'testfile2')
    error('Config must contain testfile2 to load the second file to test')
end
if ~isfield(config,'UpperLat')
    warning('Config does not contain UpperLat, no limitation of accepted Latitudes')
else
    config.UpperLat=str2double(config.UpperLat);
end
if ~isfield(config,'LowerLat')
    warning('Config does not contain LowerLat, no limitation of accepted Latitudes')
else
    config.LowerLat=str2double(config.LowerLat);
end
if ~isfield(config,'tolerance')
    warning('Config does not contain tolerance, default tolerance of 1e-8 is used')
    config.tolerance = 1e-8;
else
    config.tolerance=str2double(config.tolerance);
end



%
PAR.R_e=6378137;
PAR.f_e=1/298.257223563;
PAR.ecc_e=sqrt((2-PAR.f_e).*PAR.f_e);   % Earth Eccentricty
PAR.b_e=PAR.R_e.*sqrt(1-PAR.ecc_e.^2);  % Earth Semi Axis Minor

semimajor=6378136.3;
flattening=298.257;

facDeg2rad = pi/180;

%% Load files and find valid indices

file1 = nc_reader(config.testfile1);
file2 = nc_reader(config.testfile2);

Lat1 = file1.lat_20hz(:);
Lat2 = file2.lat_20hz(:);

% Ensure that the same points - same within a tolerance for errors - are used
tolerance = config.tolerance;
[idxLat1] = ismembertol(Lat1,Lat2,tolerance);
[idxLat2] = ismembertol(Lat2,Lat1,tolerance);

% Establish logical indices for non NaN and Values inside the Latitude
% Bounds
pos1_nan = ~isnan(Lat1);
pos2_nan = ~isnan(Lat2);

if isfield(config,'LowerLat') && isfield(config,'UpperLat')
    
    pos1_lat = (Lat1 < config.UpperLat) & (config.LowerLat < Lat1);
    pos2_lat = (Lat2 < config.UpperLat) & (config.LowerLat < Lat2);
    
    pos1 = pos1_nan & pos1_lat & idxLat1;
    pos2 = pos2_nan & pos2_lat & idxLat2;
    
else
    pos1 = pos1_nan & idxLat1;
    pos2 = pos2_nan & idxLat2;
    
end



%% Create struct from  File1

% LAT = file1.lat_20hz(:);
% pos = ~isnan(LAT);
LON = file1.lon_20hz(pos1);
LAT = Lat1(pos1);

% [x1,y1,z1]=geodetic2ecef(RADS_LAT*facDeg2rad,RADS_LON*facDeg2rad,RADS_LON*0,[semimajor, sqrt((2-(1./flattening)).*((1./flattening)))]);
% [x2,y2,z2]=geodetic2ecef(RADS_LAT*facDeg2rad,RADS_LON*facDeg2rad,RADS_LON*0,[PAR.R_e PAR.ecc_e]);

HUNC = file1.alt_cnes_20hz(pos1) - file1.range_20hz_ku(pos1);
SWH = file1.swh_20hz_ku(pos1);
U10 = get_U10(file1.sig0_20hz_ku(pos1));
sig0 = file1.sig0_20hz_ku(pos1);

structFile1 = struct('LAT', LAT,'LON', LON,...
    'pos', pos1, 'HUNC', HUNC,...
    'SWH', SWH, 'U10', U10,...
    'sig0', sig0);


%% Create struct from File2

% LAT = file2.lat_20hz(:);
% pos = ~isnan(LAT);
LON = file2.lon_20hz(pos2);
LAT = Lat2(pos2);

HUNC = file2.alt_cnes_20hz(pos2) - file2.range_20hz_ku(pos2);

SWH = file2.swh_20hz_ku(pos2);
U10 = get_U10(file2.sig0_20hz_ku(pos2));
sig0 = file2.sig0_20hz_ku(pos2);

structFile2 = struct('LAT', LAT,'LON', LON,...
    'pos', pos2, 'HUNC', HUNC,...
    'SWH', SWH, 'U10', U10,...
    'sig0', sig0);



%% Statistics
[stat,posSum,posHuncSum,structFile1,structFile2] = statisticsCalculation(structFile1,structFile2);

%% Plotting
plotting_AlongTrack(config,structFile1,structFile2,stat,posSum,posHuncSum)


end