# GRACE Crustal Deformation Analysis - Implementation Plan

## Executive Summary
Model elastic response of Earth's crust to hydrological loading using GRACE spherical harmonic coefficients (2002-2017) and validate against GPS observations from 5 California stations. This plan prioritizes existing functions and is grounded in peer-reviewed research.

## Theoretical Foundation (Peer-Reviewed)

### Core References
1. **Farrell (1972)**: "Deformation of the Earth by surface loads" - Foundational elastic loading theory
2. **Wahr et al. (1998)**: "Time variability of Earth's gravity field" - GRACE to deformation methodology
3. **Fu & Freymueller (2012)**: "Seasonal and long-term vertical deformation" - GPS-GRACE comparison
4. **Karegar et al. (2018)**: "Hybrid method for hydrologically induced vertical deformation"
5. **Wang et al. (2012)**: "Load Love numbers and Green's functions for PREM"

### Physical Basis
- **Elastic Loading Theory**: Earth responds elastically to surface mass loads on seasonal timescales (Farrell, 1972)
- **Love Numbers**: Dimensionless coefficients h_n, l_n, k_n describing elastic response (Love, 1911; Farrell, 1972)
- **Vertical Deformation Formula** (Wahr et al., 1998):
  ```
  u_r(θ,λ) = (R*ρ_w)/(3*ρ_e) * Σ_n Σ_m (h_n/(1+k_n)) * P_nm(cosθ) * [ΔC_nm*cos(mλ) + ΔS_nm*sin(mλ)]
  ```
  Where:
  - R = Earth radius (6371 km)
  - ρ_w = water density (1000 kg/m³)
  - ρ_e = average Earth density (5517 kg/m³)
  - h_n, k_n = load Love numbers for degree n

## Phase 1: Leverage Existing Functions

### 1.1 Already Available Functions
```matlab
% FROM Satellite_geodesy_assignment/:
- readSHC.m          # Reads spherical harmonic coefficients ✓
- pnm.m              # Computes normalized Legendre functions ✓
- estimate.m         # Estimates time-varying SHCs ✓

% FROM project/functions/:
- load_tenv3.m       # Loads GPS tenv3 format ✓
- fitPolynomial.m    # Polynomial detrending ✓
- solve_LSA.m        # Least squares adjustment ✓
- outlier_test.m     # Statistical outlier detection ✓
```

### 1.2 Utility Functions Available
```matlab
% Time conversion (nsea_project/):
- mjd2decyear.m      # MJD to decimal year ✓
- decyear2mjd.m      # Decimal year to MJD ✓
- mjd2date.m         # MJD to date ✓

% Statistics (lf_altimetry/):
- corr.m             # Correlation calculation ✓
- rms.m              # RMS computation ✓
- NSE.m              # Nash-Sutcliffe efficiency ✓
```

## Phase 2: New Functions Required

### 2.1 `loadLoveNumbers.m`
```matlab
function [h_n, l_n, k_n] = loadLoveNumbers(nmax, model)
% Load Love numbers based on PREM (Wang et al., 2012)
% Input:
%   nmax  - maximum degree (60 for GRACE)
%   model - 'PREM' or 'ak135' or 'iasp91'
% Output:
%   h_n, l_n, k_n - Love numbers for degrees 0:nmax
% Reference: Wang et al. (2012), Computers & Geosciences
```

### 2.2 `graceToVerticalDeformation.m`
```matlab
function [u_vertical] = graceToVerticalDeformation(cnm, snm, theta, lambda, h_n, k_n)
% Convert GRACE SH coefficients to vertical deformation
% Based on Wahr et al. (1998) equation
% Input:
%   cnm, snm - Spherical harmonic coefficients (from readSHC)
%   theta, lambda - Colatitude and longitude grids
%   h_n, k_n - Love numbers from loadLoveNumbers
% Output:
%   u_vertical - Vertical deformation in meters
% Reference: Wahr et al. (1998), JGR Solid Earth
```

### 2.3 `processGRACEfiles.m`
```matlab
function [cnm_ts, snm_ts, time_mjd] = processGRACEfiles(grace_dir, c20_file, deg1_file)
% Process all GRACE files with corrections
% Input:
%   grace_dir - Directory with .gfc files
%   c20_file - C20 replacement from SLR
%   deg1_file - Degree-1 coefficients (geocenter)
% Output:
%   cnm_ts, snm_ts - Time series of coefficients
%   time_mjd - Time stamps in MJD
% References: 
%   - Swenson et al. (2008) for degree-1
%   - Cheng et al. (2013) for C20 replacement
```

### 2.4 `extractGRACEatGPS.m`
```matlab
function [grace_at_gps] = extractGRACEatGPS(u_vertical_grid, lat_gps, lon_gps)
% Extract GRACE deformation at GPS locations
% Input:
%   u_vertical_grid - Grid from graceToVerticalDeformation
%   lat_gps, lon_gps - GPS station coordinates
% Output:
%   grace_at_gps - GRACE deformation at GPS points
% Method: Bilinear interpolation
```

### 2.5 `compareTimeSeries.m`
```matlab
function [stats] = compareTimeSeries(gps_ts, grace_ts, time_gps, time_grace)
% Statistical comparison of GPS and GRACE
% Input:
%   gps_ts, grace_ts - Time series to compare
%   time_gps, time_grace - Time vectors (MJD)
% Output:
%   stats.correlation - Pearson correlation
%   stats.rmse - Root mean square error
%   stats.nse - Nash-Sutcliffe efficiency
%   stats.amplitude_ratio - Seasonal amplitude comparison
% Reference: Fu & Freymueller (2012), JGR Solid Earth
```

## Phase 3: Main Processing Script

### 3.1 `main_grace_analysis.m`
```matlab
%% GRACE-GPS Crustal Deformation Analysis
% Based on peer-reviewed methodology from:
% - Wahr et al. (1998)
% - Fu & Freymueller (2012)
% - Karegar et al. (2018)

%% 1. Setup
clc; clear; close all;
addpath(genpath(pwd));

% Constants (IERS Conventions 2010)
R_earth = 6371000;      % m
rho_water = 1000;       % kg/m³
rho_earth = 5517;       % kg/m³
GM = 3.986004418e14;    % m³/s²

%% 2. Load GPS Data
stations = {'P056', 'P140', 'P245', 'P308', 'P345'};
gps_data = struct();
for i = 1:length(stations)
    gps_data.(stations{i}) = load_tenv3([stations{i} '.tenv3']);
    % Detrend using existing fitPolynomial.m
    [gps_detrended] = fitPolynomial(gps_data.(stations{i}).t, ...
                                   gps_data.(stations{i}).up, 2, ...
                                   mean(gps_data.(stations{i}).sig_um));
end

%% 3. Load and Process GRACE
[cnm_ts, snm_ts, time_mjd] = processGRACEfiles('data/grace/', ...
                                                'data/aux/C20_RL05.txt', ...
                                                'data/aux/deg1_coef.txt');

%% 4. Load Love Numbers (PREM)
nmax = 60;  % Maximum degree from GRACE files
[h_n, l_n, k_n] = loadLoveNumbers(nmax, 'PREM');

%% 5. Convert GRACE to Deformation
% Create spatial grid
res = 0.5;  % degrees
lat = -90:res:90;
lon = -180:res:180;
[LON, LAT] = meshgrid(lon, lat);
theta = (90 - LAT) * pi/180;  % colatitude

% Process each month
n_months = size(cnm_ts, 3);
u_vertical_ts = zeros(length(lat), length(lon), n_months);

for month = 1:n_months
    % Use existing pnm.m for Legendre functions
    u_vertical_ts(:,:,month) = graceToVerticalDeformation(...
        cnm_ts(:,:,month), snm_ts(:,:,month), theta, LON*pi/180, h_n, k_n);
end

%% 6. Extract at GPS Locations and Compare
gps_coords = load('data/gps/GPSLatLong.tenv3');
results = struct();

for i = 1:length(stations)
    % Extract GRACE at GPS location
    grace_at_gps = extractGRACEatGPS(u_vertical_ts, ...
                                     gps_coords.lat(i), ...
                                     gps_coords.lon(i));
    
    % Compare time series
    results.(stations{i}) = compareTimeSeries(...
        gps_detrended.(stations{i}).v_detrended, ...
        grace_at_gps, ...
        gps_data.(stations{i}).t, ...
        time_mjd);
end

%% 7. Generate Outputs
generateReport(results, stations);
saveas(gcf, 'output/grace_gps_comparison.png');
save('output/results.mat', 'results', 'u_vertical_ts', 'time_mjd');
```

## Phase 4: Validation & Quality Control

### 4.1 Physical Validation
- **Amplitude Check**: Seasonal deformation should be 1-10 mm (Fu & Freymueller, 2012)
- **Phase Check**: Maximum loading typically in winter/spring (Argus et al., 2014)
- **Spatial Coherence**: Nearby stations should show similar patterns

### 4.2 Statistical Validation
- **Correlation Threshold**: R > 0.5 for seasonal signals (typical in literature)
- **RMSE Target**: < 5 mm for vertical component
- **NSE Target**: > 0.3 (acceptable model performance)

### 4.3 Known Issues & Solutions
- **Geocenter Motion**: Add degree-1 terms (Swenson et al., 2008)
- **C20 Accuracy**: Replace with SLR values (Cheng et al., 2013)
- **Spatial Resolution**: GRACE ~300 km vs GPS point → expect smoothing
- **Glacial Isostatic Adjustment**: May need GIA correction (Peltier et al., 2018)

## Phase 5: Expected Scientific Outcomes

### 5.1 Primary Results
1. **Correlation Maps**: GPS vs GRACE vertical deformation correlation
2. **Seasonal Analysis**: Amplitude and phase of hydrological loading
3. **Time Series**: Monthly deformation 2002-2017
4. **Statistical Metrics**: RMSE, correlation, NSE for each station

### 5.2 Physical Interpretation
- **Hydrological Loading**: Dominant seasonal signal in California
- **Drought Impact**: Reduced loading during 2012-2016 California drought
- **Spatial Patterns**: Sierra Nevada snowpack influence

## Phase 6: Computational Considerations

### 6.1 Memory Management
```matlab
% Process GRACE files in batches to avoid memory issues
batch_size = 12;  % months
clear large_arrays after use
```

### 6.2 Performance Optimization
- Use existing `pnm.m` for Legendre functions (already optimized)
- Vectorize operations where possible
- Pre-allocate arrays

## Timeline
1. **Week 1**: Implement new functions (2.1-2.5)
2. **Week 2**: Process GRACE data with corrections
3. **Week 3**: GPS data processing and detrending
4. **Week 4**: Comparison and validation
5. **Week 5**: Generate figures and report

## References

### Primary Literature
1. Farrell, W.E. (1972). Deformation of the Earth by surface loads. Reviews of Geophysics, 10(3), 761-797.
2. Wahr, J., Molenaar, M., & Bryan, F. (1998). Time variability of the Earth's gravity field. JGR Solid Earth, 103(B12), 30205-30229.
3. Fu, Y., & Freymueller, J.T. (2012). Seasonal and long‐term vertical deformation in the Nepal Himalaya. JGR Solid Earth, 117(B3).
4. Karegar, M.A., et al. (2018). A new hybrid method for estimating hydrologically induced vertical deformation. Journal of Advances in Modeling Earth Systems, 10(5), 1196-1217.

### Supporting Literature
5. Wang, H., et al. (2012). Load Love numbers and Green's functions for elastic Earth models. Computers & Geosciences, 49, 190-199.
6. Swenson, S., et al. (2008). Estimating geocenter variations from GRACE and ocean model. JGR Solid Earth, 113(B8).
7. Cheng, M., et al. (2013). Variations in the Earth's oblateness during the past 28 years. JGR Solid Earth, 118(2), 740-747.
8. Dziewonski, A.M., & Anderson, D.L. (1981). Preliminary reference Earth model. Physics of the Earth and Planetary Interiors, 25(4), 297-356.