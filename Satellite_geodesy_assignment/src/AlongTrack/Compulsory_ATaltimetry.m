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
addpath(genpath('../../lib/matlab/'));

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

%% Task 2 interpolate at 1 Hz

% Load files
file1 = nc_reader(config.testfile1);
file2 = nc_reader(config.testfile2);

% Extract 20 Hz data from file1
Lat1 = file1.lat_20hz(:);
Lon1 = file1.lon_20hz(:);
HUNC1 = file1.alt_cnes_20hz(:) - file1.range_20hz_ku(:);
SWH1 = file1.swh_20hz_ku(:);
U10_1 = get_U10(file1.sig0_20hz_ku(:));
sig0_1 = file1.sig0_20hz_ku(:);

% Extract 20 Hz data from file2
Lat2 = file2.lat_20hz(:);
Lon2 = file2.lon_20hz(:);
HUNC2 = file2.alt_cnes_20hz(:) - file2.range_20hz_ku(:);
SWH2 = file2.swh_20hz_ku(:);
U10_2 = get_U10(file2.sig0_20hz_ku(:));
sig0_2 = file2.sig0_20hz_ku(:);

% Generate time vectors (assuming uniform time steps)
t_20hz_1 = (0:length(HUNC1)-1) / 20; % 20 Hz time axis
t_20hz_2 = (0:length(HUNC2)-1) / 20; % 20 Hz time axis
t_1hz = (0:max(t_20hz_1))'; % 1 Hz time axis

% Interpolate data to 1 Hz for file1
HUNC1_1hz = interp1(t_20hz_1, HUNC1, t_1hz, 'linear');
SWH1_1hz = interp1(t_20hz_1, SWH1, t_1hz, 'linear');
U10_1hz = interp1(t_20hz_1, U10_1, t_1hz, 'linear');
sig0_1hz = interp1(t_20hz_1, sig0_1, t_1hz, 'linear');
Lat1_1hz = interp1(t_20hz_1, Lat1, t_1hz, 'linear');
Lon1_1hz = interp1(t_20hz_1, Lon1, t_1hz, 'linear');

% Interpolate data to 1 Hz for file2
HUNC2_1hz = interp1(t_20hz_2, HUNC2, t_1hz, 'linear');
SWH2_1hz = interp1(t_20hz_2, SWH2, t_1hz, 'linear');
U10_2hz = interp1(t_20hz_2, U10_2, t_1hz, 'linear');
sig0_2hz = interp1(t_20hz_2, sig0_2, t_1hz, 'linear');
Lat2_1hz = interp1(t_20hz_2, Lat2, t_1hz, 'linear');
Lon2_1hz = interp1(t_20hz_2, Lon2, t_1hz, 'linear');


structFile1_1hz = struct('LAT', Lat1_1hz, 'LON', Lon1_1hz, ...
    'HUNC', HUNC1_1hz, 'SWH', SWH1_1hz, 'U10', U10_1hz, 'sig0', sig0_1hz);

structFile2_1hz = struct('LAT', Lat2_1hz, 'LON', Lon2_1hz, ...
    'HUNC', HUNC2_1hz, 'SWH', SWH2_1hz, 'U10', U10_2hz, 'sig0', sig0_2hz);

[stat_1hz,posSum_1hz,posHuncSum_1hz,structFile1_1hz,structFile2_1hz] = statisticsCalculation(structFile1_1hz,structFile2_1hz);
config_1hz = config;
config_1hz.nameFile1 = 'SINC2_1hz';
config_1hz.nameFile2 = 'TALES_1hz';
% Plotting
plotting_AlongTrack(config_1hz,structFile1_1hz,structFile2_1hz,stat_1hz,posSum_1hz,posHuncSum_1hz)



%% Task 3 Add corrections to study sea level
% Extract 1 Hz corrections from file1
dry_tropo1_1hz = file1.dry_tropo(:);
wet_tropo1_1hz = file1.wet_tropo(:);
iono1_1hz = file1.iono(:);  
ssb1_1hz = file1.ssb(:);
inv_bar1_1hz = file1.inv_bar(:);  % Inverse Barometer Correction

% Tide corrections (1 Hz)
tide_solid1_1hz = file1.tide_solid(:);
tide_ocean1_1hz = file1.tide_ocean(:);
tide_load1_1hz = file1.tide_load(:);
tide_pole1_1hz = file1.tide_pole(:);
tide_equil1_1hz = file1.tide_equil(:);
tide_non_equil1_1hz = file1.tide_non_equil(:);

% Compute total correction for file1 (already in 1 Hz)
C_total1_1hz = dry_tropo1_1hz + wet_tropo1_1hz + iono1_1hz +inv_bar1_1hz+... %ssb1_1hz + ...
               (tide_solid1_1hz + tide_ocean1_1hz + tide_load1_1hz + tide_pole1_1hz + ...
               tide_equil1_1hz + tide_non_equil1_1hz);

% Compute corrected SSH for file1
SSH_corrected1_1hz = HUNC1_1hz - C_total1_1hz; 

% Repeat for file2
dry_tropo2_1hz = file2.dry_tropo(:);
wet_tropo2_1hz = file2.wet_tropo(:);
iono2_1hz = file2.iono(:);
ssb2_1hz = file2.ssb(:);
inv_bar2_1hz = file2.inv_bar(:);  % Inverse Barometer Correction

tide_solid2_1hz = file2.tide_solid(:);
tide_ocean2_1hz = file2.tide_ocean(:);
tide_load2_1hz = file2.tide_load(:);
tide_pole2_1hz = file2.tide_pole(:);
tide_equil2_1hz = file2.tide_equil(:);
tide_non_equil2_1hz = file2.tide_non_equil(:);

C_total2_1hz = dry_tropo2_1hz + wet_tropo2_1hz + iono2_1hz +inv_bar2_1hz+...% ssb2_1hz + ...
               (tide_solid2_1hz + tide_ocean2_1hz + tide_load2_1hz + tide_pole2_1hz + ...
               tide_equil2_1hz + tide_non_equil2_1hz);

SSH_corrected2_1hz = HUNC2_1hz - C_total2_1hz;

structFile1_1hz_c = struct('LAT', Lat1_1hz, 'LON', Lon1_1hz, ...
    'HUNC', SSH_corrected1_1hz, 'SWH', SWH1_1hz, 'U10', U10_1hz, 'sig0', sig0_1hz);

structFile2_1hz_c = struct('LAT', Lat2_1hz, 'LON', Lon2_1hz, ...
    'HUNC', SSH_corrected2_1hz, 'SWH', SWH2_1hz, 'U10', U10_2hz, 'sig0', sig0_2hz);

config_1hz_c = config;
config_1hz_c.nameFile1 = 'SINC2_1hz_c';
config_1hz_c.nameFile2 = 'TALES_1hz_c';

[stat_1hz_c,posSum_1hz_c,posHuncSum_1hz_c,structFile1_1hz_c,structFile2_1hz_c] = statisticsCalculation(structFile1_1hz_c,structFile2_1hz_c);

% Plotting
plotting_AlongTrack(config_1hz_c,structFile1_1hz_c,structFile2_1hz_c,stat_1hz_c,posSum_1hz_c,posHuncSum_1hz_c)


%% Task 4 Add corrections to compare with tide gauge

% Extract additional 1 Hz corrections from file1
geoid1_20hz = file1.geoid(:);  
% Interpolate geoid to 1 Hz
geoid1_1hz = interp1(t_20hz_1, geoid1_20hz, t_1hz, 'linear');

% Compute geoid-corrected SSH
SSH_geoid_corrected1_1hz = structFile1_1hz_c.HUNC - geoid1_1hz;


% Store in a new struct for tide gauge validation
structFile1_1hz_gauge = struct('LAT', structFile1_1hz_c.LAT, 'LON', structFile1_1hz_c.LON, ...
    'HUNC', SSH_geoid_corrected1_1hz, 'SWH', structFile1_1hz_c.SWH, ...
    'U10', structFile1_1hz_c.U10, 'sig0', structFile1_1hz_c.sig0);

% Repeat for file2
geoid2_20hz = file2.geoid(:);  
geoid2_1hz = interp1(t_20hz_2, geoid2_20hz, t_1hz, 'linear');

SSH_geoid_corrected2_1hz = structFile2_1hz_c.HUNC - geoid2_1hz;

structFile2_1hz_gauge = struct('LAT', structFile2_1hz_c.LAT, 'LON', structFile2_1hz_c.LON, ...
    'HUNC', SSH_geoid_corrected2_1hz, 'SWH', structFile2_1hz_c.SWH, ...
    'U10', structFile2_1hz_c.U10, 'sig0', structFile2_1hz_c.sig0);

config_1hz_gauge = config;
config_1hz_gauge.nameFile1 = 'SINC2_1hz_gauge';
config_1hz_gauge.nameFile2 = 'TALES_1hz_gauge';

[stat_1hz_gauge,posSum_1hz_gauge,posHuncSum_1hz_gauge,structFile1_1hz_gauge,structFile2_1hz_gauge] = statisticsCalculation(structFile1_1hz_gauge,structFile2_1hz_gauge);

% Plotting
plotting_AlongTrack(config_1hz_gauge,structFile1_1hz_gauge,structFile2_1hz_gauge,stat_1hz_gauge,posSum_1hz_gauge,posHuncSum_1hz_gauge)

%% Task 5 
base_directory = '/Users/jaradata/Documents/Hamza_MSc/Satellite_geodesy_assignment/data/Compulsory_timeseries';
[dahiti_combined, theia_combined] = read_and_combine_time_series(base_directory);

dahiti_cleaned = rmmissing(dahiti_combined);
theia_cleaned = rmmissing(theia_combined);
    
% Sort the data by date
dahiti_cleaned = sortrows(dahiti_cleaned, 'Date');
theia_cleaned = sortrows(theia_cleaned, 'Date');

figure;
hold on; 
% Plot Dahiti data
errorbar(dahiti_cleaned.Date, dahiti_cleaned.WaterLevel, dahiti_cleaned.Error, 'bo', 'DisplayName', 'Dahiti', 'MarkerSize', 12, 'LineWidth', 3);
% Plot Theia data
errorbar(theia_cleaned.Date, theia_cleaned.WaterLevel, theia_cleaned.Error, 'ro', 'DisplayName', 'Theia', 'MarkerSize', 12, 'LineWidth', 3);
% Add labels and legend
xlabel('Date');
ylabel('Water Level (m)');
title('Comparison of Dahiti and Theia Time Series');
legend('show');
 ax=gca;ax.XAxis.FontSize = 22;ax.YAxis.FontSize = 22;ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;ax.Title.FontSize = 24;ax.FontSize=18;
grid on;
hold off;
end