# Part 2: GPS Time Series Analysis
## Extracting Vertical Crustal Deformation from GPS Observations

---

## Table of Contents
1. [Overview and Principles](#1-overview-and-principles)
2. [GPS Theory and Positioning](#2-gps-theory-and-positioning)
3. [TENV3 Data Format](#3-tenv3-data-format)
4. [Time Series Analysis Theory](#4-time-series-analysis-theory)
5. [Step-by-Step Processing Pipeline](#5-step-by-step-processing-pipeline)
6. [Statistical Methods](#6-statistical-methods)
7. [Quality Control and Validation](#7-quality-control-and-validation)
8. [Error Analysis](#8-error-analysis)

---

## 1. Overview and Principles

### 1.1 GPS for Crustal Deformation

The Global Positioning System (GPS) provides precise measurements of surface positions that can detect crustal deformation with millimeter-level accuracy. When combined with GRACE satellite gravity data, GPS observations serve as an independent validation technique for elastic loading theory.

**Key Advantages of GPS**:
- **High Temporal Resolution**: Daily or sub-daily measurements
- **Point Precision**: Millimeter-level accuracy at individual stations
- **Direct Measurement**: Measures actual surface displacements (not inferential)
- **Long Time Series**: Multi-decade records available at many stations

**Applications to Hydrological Loading**:
GPS stations detect elastic crustal deformation caused by:
- Seasonal water storage variations
- Snow and ice loading/unloading
- Soil moisture changes
- Groundwater variations

### 1.2 Physical Basis

**Elastic Response Principle**:
When surface loads change (e.g., seasonal hydrological variations), the solid Earth responds elastically, causing:
- **Vertical Displacements**: Typically 1-10 mm for hydrological loading
- **Horizontal Displacements**: Generally smaller than vertical (ratio ~0.6)
- **Tilt Changes**: Angular deformation measured by GPS receiver positioning

**Loading vs. GPS Response**:
```
Surface Load Increase → Crustal Subsidence (negative vertical)
Surface Load Decrease → Crustal Uplift (positive vertical)
```

### 1.3 Integration with GRACE Analysis

GPS time series analysis (Part 2) provides the observational foundation for comparing with GRACE-based theoretical predictions (Part 1). This enables validation of:
- Elastic loading theory accuracy
- Love number models (PREM)
- GRACE data quality and processing

---

## 2. GPS Theory and Positioning

### 2.1 GPS Positioning Fundamentals

**Basic GPS Principle**:
GPS positioning determines receiver location by measuring distances to multiple satellites with known positions.

**Distance Measurement**:
```
ρ = c × Δt
```
Where:
- `ρ` = pseudorange to satellite [m]
- `c` = speed of light = 2.998 × 10⁸ m/s
- `Δt` = travel time difference [s]

**Position Solution**:
For position `(x, y, z)` and receiver clock bias `b`, solve the system:

```
ρᵢ = √[(x-xᵢ)² + (y-yᵢ)² + (z-zᵢ)²] + c·b + εᵢ
```

Where:
- `(xᵢ, yᵢ, zᵢ)` = known satellite i position
- `εᵢ` = measurement errors and biases

### 2.2 Precision Positioning Techniques

**Differential GPS (DGPS)**:
Uses reference stations to eliminate common errors:
- Atmospheric delays
- Satellite clock errors  
- Orbital errors

**Real-Time Kinematic (RTK)**:
Carrier phase-based positioning achieving centimeter accuracy:

```
φ = (1/λ) × ρ + N + ε
```

Where:
- `φ` = carrier phase measurement [cycles]
- `λ` = carrier wavelength [m]
- `N` = integer ambiguity [cycles]
- `ε` = phase errors [cycles]

**Precise Point Positioning (PPP)**:
Single-receiver technique using precise satellite products:
- Precise orbits and clocks
- Ionospheric and tropospheric corrections
- Achieves millimeter accuracy in post-processing

### 2.3 Coordinate Systems and Reference Frames

**World Geodetic System 1984 (WGS84)**:
- Geocentric coordinate system
- Used by GPS satellites
- Earth-Centered, Earth-Fixed (ECEF) frame

**International Terrestrial Reference Frame (ITRF)**:
- Scientific reference frame for precise applications
- Accounts for tectonic plate motion
- Updated periodically (ITRF2014, ITRF2020)

**Coordinate Transformations**:
```matlab
% ECEF to Geographic (WGS84)
[lat, lon, h] = ecef2geodetic(X, Y, Z, wgs84Ellipsoid);

% Geographic to Local Tangent Plane
[east, north, up] = geodetic2enu(lat, lon, h, lat0, lon0, h0, wgs84Ellipsoid);
```

### 2.4 Error Sources in GPS Positioning

**Satellite-Related Errors**:
- Orbital errors: ~2-5 m (corrected with precise ephemerides)
- Clock errors: ~2 m equivalent (corrected with precise clocks)

**Atmospheric Errors**:
- Ionospheric delay: 1-50 m (dual-frequency correction)
- Tropospheric delay: 2-25 m (modeling and estimation)

**Receiver and Site Errors**:
- Multipath: mm to cm level
- Antenna phase center variations: mm level
- Monument stability: mm to cm level

---

## 3. TENV3 Data Format

### 3.1 Format Specification

The TENV3 format is widely used for GPS time series data, particularly from the Nevada Geodetic Laboratory (NGL) and similar processing centers.

**File Structure**:
```
# Header information
# Site: P140  Lat:   37.2345  Lon: -115.1234  Height:   1234.567
# Reference position and velocity information
# Column definitions...

SITE YR MO DA HR MN SS DECYEAR      MJD       UP         E          N        UU        UN        UE        EU        NN        EN       UUV       UNV       UEV       NU        EU        UU
P140 02  4  5 12  0  0 2002.2603 52369.5000  0.0045  -0.0012   0.0023  0.0015  0.0008  0.0010  0.0008  0.0012  0.0009   1.5e-06  0.8e-06  1.0e-06  0.8e-06  1.2e-06  1.5e-06
P140 02  4  6 12  0  0 2002.2630 52370.5000  0.0047  -0.0010   0.0025  0.0015  0.0008  0.0010  0.0008  0.0012  0.0009   1.5e-06  0.8e-06  1.0e-06  0.8e-06  1.2e-06  1.5e-06
...
```

**Column Definitions**:
1. **SITE**: 4-character station identifier
2. **YR**: Year (2-digit)  
3. **MO**: Month (1-12)
4. **DA**: Day of month (1-31)
5. **HR**: Hour (0-23)
6. **MN**: Minute (0-59)
7. **SS**: Second (0-59)
8. **DECYEAR**: Decimal year
9. **MJD**: Modified Julian Date
10. **UP**: Vertical component [m]
11. **E**: East component [m]
12. **N**: North component [m]
13. **UU**: Up-Up covariance [m²]
14. **UN**: Up-North covariance [m²]
15. **UE**: Up-East covariance [m²]
16. **EU**: East-Up covariance [m²]  
17. **NN**: North-North covariance [m²]
18. **EN**: East-North covariance [m²]

### 3.2 Data Loading Implementation

**MATLAB Loading Function**:
```matlab
function gps_data = load_tenv3(filename)
    % Load TENV3 format GPS time series file
    
    % Check file existence
    if ~exist(filename, 'file')
        error('File not found: %s', filename);
    end
    
    % Define format specification
    formatSpec = '%4s%3d%3d%3d%3d%3d%3d%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%10f%15f%15f%15f%15f%15f%15f%*[^\n]';
    
    % Open file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    % Read header (skip comment lines starting with #)
    header_lines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, '#')
            header_lines{end+1} = line;
        else
            % Reset file pointer to beginning of data
            frewind(fid);
            break;
        end
    end
    
    % Read data
    try
        dataArray = textscan(fid, formatSpec, 'CommentStyle', '#', ...
                           'EmptyValue', NaN, 'ReturnOnError', false);
        fclose(fid);
    catch ME
        fclose(fid);
        error('Error reading TENV3 file: %s', ME.message);
    end
    
    % Extract data into structure
    gps_data = struct();
    
    % Time information
    gps_data.site = dataArray{1};
    gps_data.year = dataArray{2} + 2000; % Convert 2-digit to 4-digit year
    gps_data.month = dataArray{3};
    gps_data.day = dataArray{4};
    gps_data.decyear = dataArray{8};
    gps_data.mjd = dataArray{9};
    
    % Position components
    gps_data.up = dataArray{10};
    gps_data.east = dataArray{11}; 
    gps_data.north = dataArray{12};
    
    % Uncertainties (diagonal terms)
    gps_data.sig_up = sqrt(abs(dataArray{13}));    % Standard deviation from variance
    gps_data.sig_east = sqrt(abs(dataArray{17}));
    gps_data.sig_north = sqrt(abs(dataArray{18}));
    
    % Covariance terms (for full error propagation if needed)
    gps_data.cov_up_up = dataArray{13};
    gps_data.cov_up_north = dataArray{14};
    gps_data.cov_up_east = dataArray{15};
    gps_data.cov_east_up = dataArray{16};
    gps_data.cov_north_north = dataArray{17};
    gps_data.cov_east_north = dataArray{18};
    
    % Data quality check
    n_obs = length(gps_data.mjd);
    fprintf('Loaded %d observations from %s\n', n_obs, filename);
    
    % Remove invalid observations
    valid_idx = isfinite(gps_data.up) & isfinite(gps_data.east) & isfinite(gps_data.north) & ...
                isfinite(gps_data.mjd) & gps_data.mjd > 0;
    
    if sum(valid_idx) < n_obs
        fprintf('Removed %d invalid observations\n', n_obs - sum(valid_idx));
        
        % Filter all arrays
        fields = fieldnames(gps_data);
        for i = 1:length(fields)
            if length(gps_data.(fields{i})) == n_obs
                gps_data.(fields{i}) = gps_data.(fields{i})(valid_idx);
            end
        end
    end
    
    % Summary statistics
    fprintf('Time span: %.3f - %.3f (%.1f years)\n', ...
            min(gps_data.decyear), max(gps_data.decyear), ...
            max(gps_data.decyear) - min(gps_data.decyear));
    
    fprintf('Vertical component statistics:\n');
    fprintf('  Mean: %.3f mm\n', mean(gps_data.up) * 1000);
    fprintf('  Std:  %.3f mm\n', std(gps_data.up) * 1000);
    fprintf('  Range: %.3f to %.3f mm\n', ...
            min(gps_data.up) * 1000, max(gps_data.up) * 1000);
end
```

### 3.3 Coordinate File Format

**GPS Station Coordinates (GPSLatLong.tenv3)**:
```
# Station coordinates in WGS84
# SITE      LAT        LON       HEIGHT
P056   34.1234   -118.5678   123.456
P140   37.2345   -115.1234   1234.567
P245   39.3456   -119.8765   1456.789
...
```

**Loading Station Coordinates**:
```matlab
function [station_info] = loadGPScoordinates(coord_file)
    % Load GPS station coordinate information
    
    fid = fopen(coord_file, 'r');
    if fid == -1
        error('Cannot open coordinate file: %s', coord_file);
    end
    
    % Read data (skip header lines)
    data = textscan(fid, '%s %f %f %f', 'CommentStyle', '#');
    fclose(fid);
    
    % Extract information
    station_info = struct();
    station_info.sites = data{1};
    station_info.latitude = data{2};
    station_info.longitude = data{3};
    station_info.height = data{4};
    station_info.n_stations = length(station_info.sites);
    
    fprintf('Loaded coordinates for %d GPS stations\n', station_info.n_stations);
    
    % Validation
    assert(all(abs(station_info.latitude) <= 90), 'Invalid latitude values');
    assert(all(abs(station_info.longitude) <= 180), 'Invalid longitude values');
end
```

---

## 4. Time Series Analysis Theory

### 4.1 Mathematical Model for GPS Time Series

GPS time series contain multiple signal components that must be separated to extract the hydrological loading signal:

```
y(t) = a₀ + a₁t + a₂t² + Σᵢ[Aᵢcos(2πfᵢt) + Bᵢsin(2πfᵢt)] + η(t) + ε(t)
```

Where:
- `y(t)` = observed GPS position at time t
- `a₀ + a₁t + a₂t²` = polynomial trend (secular motion + acceleration)
- `Aᵢ, Bᵢ` = amplitudes of periodic signals at frequency fᵢ
- `η(t)` = hydrological loading signal (target of interest)
- `ε(t)` = observation noise and unmodeled effects

### 4.2 Signal Components

**4.2.1 Secular Trends**

**Linear Trend**: Represents constant velocity (tectonic motion):
```
v = a₁ = dy/dt [m/year]
```

**Quadratic Trend**: Represents acceleration (postseismic relaxation):
```
a = 2a₂ [m/year²]
```

**Physical Sources**:
- Tectonic plate motion
- Glacial isostatic adjustment (GIA)
- Postseismic relaxation after earthquakes
- Long-term anthropogenic effects

**4.2.2 Periodic Signals**

**Annual Signal** (f = 1.0 cycle/year):
```
y_annual(t) = A_annual·cos(2πt) + B_annual·sin(2πt)
```

**Semi-annual Signal** (f = 2.0 cycles/year):
```  
y_semi(t) = A_semi·cos(4πt) + B_semi·sin(4πt)
```

**Amplitude and Phase Representation**:
```
A·cos(2πft) + B·sin(2πft) = R·cos(2πft - φ)
```

Where:
- `R = √(A² + B²)` = amplitude
- `φ = arctan(B/A)` = phase [radians]

**Physical Sources of Periodic Signals**:
- **Annual**: Primarily hydrological loading (target signal)
- **Semi-annual**: Mixed hydrological and atmospheric effects
- **Daily/Sub-daily**: Atmospheric and ocean tides (usually filtered)

### 4.3 Detrending Methods

**4.3.1 Polynomial Detrending**

**Linear Detrending (p=1)**:
```matlab
function [detrended_data, trend_params] = polynomialDetrend(time, data, degree, weights)
    % Fit polynomial trend using weighted least squares
    
    % Create design matrix
    A = zeros(length(time), degree+1);
    for p = 0:degree
        A(:, p+1) = (time - time(1)).^p; % Center time at start
    end
    
    % Weight matrix
    if nargin < 4 || isempty(weights)
        W = eye(length(time));
    else
        W = diag(weights);
    end
    
    % Weighted least squares solution
    trend_params = (A' * W * A) \ (A' * W * data);
    
    % Compute trend
    trend = A * trend_params;
    
    % Detrended data
    detrended_data = data - trend;
end
```

**Quadratic Detrending (p=2)**:
Includes acceleration term for postseismic analysis:

```matlab
% Design matrix for quadratic trend
A = [ones(n,1), (time-time(1)), (time-time(1)).^2];
params = A \ data; % Least squares solution

% Trend components
a0 = params(1); % Offset
a1 = params(2); % Velocity [m/year]  
a2 = params(3); % Acceleration [m/year²]
```

**4.3.2 Robust Detrending**

For data with outliers, use robust regression methods:

```matlab
function [detrended_data, trend_params] = robustDetrend(time, data, degree)
    % Robust polynomial detrending using iteratively reweighted least squares
    
    % Create design matrix
    A = zeros(length(time), degree+1);
    for p = 0:degree
        A(:, p+1) = (time - time(1)).^p;
    end
    
    % Initial least squares solution
    params = A \ data;
    
    % Iteratively reweight
    for iter = 1:5
        % Compute residuals
        residuals = data - A * params;
        
        % Robust weights (Huber weights)
        mad_residual = median(abs(residuals - median(residuals)));
        normalized_residuals = residuals / (1.4826 * mad_residual);
        
        weights = ones(size(residuals));
        large_residual = abs(normalized_residuals) > 1.345;
        weights(large_residual) = 1.345 ./ abs(normalized_residuals(large_residual));
        
        % Weighted least squares update
        W = diag(weights);
        params = (A' * W * A) \ (A' * W * data);
    end
    
    % Final results
    trend_params = params;
    detrended_data = data - A * params;
end
```

### 4.4 Harmonic Analysis

**4.4.1 Least Squares Harmonic Fitting**

Extract periodic components using least squares:

```matlab
function [harmonic_params, harmonic_signal] = fitHarmonics(time, data, frequencies)
    % Fit harmonic components at specified frequencies
    
    n_freq = length(frequencies);
    n_obs = length(time);
    
    % Create design matrix
    A = zeros(n_obs, 2*n_freq + 1);
    A(:, 1) = ones(n_obs, 1); % Constant term
    
    for i = 1:n_freq
        f = frequencies(i);
        A(:, 2*i) = cos(2*pi*f*time);     % Cosine term
        A(:, 2*i+1) = sin(2*pi*f*time);   % Sine term
    end
    
    % Least squares solution
    harmonic_params = A \ data;
    
    % Reconstructed harmonic signal
    harmonic_signal = A * harmonic_params;
    
    % Extract amplitude and phase for each frequency
    amplitude = zeros(n_freq, 1);
    phase = zeros(n_freq, 1);
    
    for i = 1:n_freq
        cos_coeff = harmonic_params(2*i);
        sin_coeff = harmonic_params(2*i+1);
        
        amplitude(i) = sqrt(cos_coeff^2 + sin_coeff^2);
        phase(i) = atan2(sin_coeff, cos_coeff);
    end
    
    % Store results in structure
    results = struct();
    results.constant = harmonic_params(1);
    results.frequencies = frequencies;
    results.cos_coeffs = harmonic_params(2:2:end);
    results.sin_coeffs = harmonic_params(3:2:end);
    results.amplitudes = amplitude;
    results.phases = phase;
    
    harmonic_params = results;
end
```

**4.4.2 Spectral Analysis**

Use Fourier analysis to identify dominant frequencies:

```matlab
function [frequencies, power_spectrum] = computePowerSpectrum(time, data)
    % Compute power spectrum to identify dominant frequencies
    
    % Ensure uniform time sampling (interpolate if necessary)
    dt = median(diff(time));
    time_uniform = time(1):dt:time(end);
    data_uniform = interp1(time, data, time_uniform, 'linear');
    
    % Remove NaN values
    valid_idx = ~isnan(data_uniform);
    time_uniform = time_uniform(valid_idx);
    data_uniform = data_uniform(valid_idx);
    
    % Apply window to reduce spectral leakage
    window = hanning(length(data_uniform));
    data_windowed = data_uniform .* window';
    
    % Compute FFT
    Y = fft(data_windowed);
    power = abs(Y).^2 / length(Y);
    
    % Create frequency vector
    fs = 1/dt; % Sampling frequency [1/years]
    frequencies = (0:length(Y)-1) * fs / length(Y);
    
    % Keep only positive frequencies
    half_length = floor(length(Y)/2) + 1;
    frequencies = frequencies(1:half_length);
    power_spectrum = power(1:half_length);
    
    % Scale power (account for two-sided spectrum)
    power_spectrum(2:end-1) = 2 * power_spectrum(2:end-1);
end
```

---

## 5. Step-by-Step Processing Pipeline

### 5.1 Pipeline Overview

The GPS processing pipeline converts raw TENV3 time series into detrended vertical displacement records suitable for comparison with GRACE:

```
Step 1: Load GPS Time Series Data
Step 2: Quality Control and Outlier Detection
Step 3: Gap Analysis and Data Continuity Check
Step 4: Polynomial Detrending
Step 5: Harmonic Analysis (Optional)
Step 6: Reference Frame Corrections
Step 7: Uncertainty Propagation
Step 8: Output Formatting for GRACE Comparison
```

### 5.2 Step 1: Load GPS Time Series Data

**Objective**: Load GPS station coordinates and time series data from TENV3 files.

**Implementation**:
```matlab
function gps_database = loadAllGPSstations(gps_data_dir, coord_file)
    % Load complete GPS database for analysis
    
    % Load station coordinates
    station_info = loadGPScoordinates(coord_file);
    n_stations = station_info.n_stations;
    
    % Initialize storage
    gps_database = struct();
    gps_database.station_info = station_info;
    gps_database.time_series = cell(n_stations, 1);
    gps_database.valid_stations = false(n_stations, 1);
    
    fprintf('Loading time series for %d GPS stations...\n', n_stations);
    
    % Process each station
    for i = 1:n_stations
        site_name = station_info.sites{i};
        
        % Construct filename (try multiple naming conventions)
        possible_files = {
            fullfile(gps_data_dir, sprintf('%s.tenv3', site_name)),
            fullfile(gps_data_dir, sprintf('%s.tenv', site_name)),
            fullfile(gps_data_dir, sprintf('%s_clean.tenv3', site_name))
        };
        
        time_series = [];
        for j = 1:length(possible_files)
            if exist(possible_files{j}, 'file')
                try
                    time_series = load_tenv3(possible_files{j});
                    break;
                catch ME
                    fprintf('Error loading %s: %s\n', possible_files{j}, ME.message);
                end
            end
        end
        
        if ~isempty(time_series)
            gps_database.time_series{i} = time_series;
            gps_database.valid_stations(i) = true;
            
            if mod(i, 10) == 0 || i == n_stations
                fprintf('  Loaded %d/%d stations\n', i, n_stations);
            end
        else
            fprintf('Warning: No data file found for station %s\n', site_name);
        end
    end
    
    n_valid = sum(gps_database.valid_stations);
    fprintf('Successfully loaded %d/%d GPS stations\n', n_valid, n_stations);
    
    % Filter to valid stations only
    gps_database.station_info.sites = station_info.sites(gps_database.valid_stations);
    gps_database.station_info.latitude = station_info.latitude(gps_database.valid_stations);
    gps_database.station_info.longitude = station_info.longitude(gps_database.valid_stations);
    gps_database.station_info.height = station_info.height(gps_database.valid_stations);
    gps_database.time_series = gps_database.time_series(gps_database.valid_stations);
    gps_database.station_info.n_stations = n_valid;
end
```

### 5.3 Step 2: Quality Control and Outlier Detection

**Objective**: Identify and remove outliers, bad observations, and unrealistic values.

**Statistical Outlier Detection**:
```matlab
function [data_clean, outlier_flags] = detectOutliers(time, data, method, threshold)
    % Detect outliers using various statistical methods
    
    if nargin < 3, method = 'iqr'; end
    if nargin < 4, threshold = 3; end
    
    n_obs = length(data);
    outlier_flags = false(n_obs, 1);
    
    switch method
        case 'iqr'
            % Interquartile Range method
            Q1 = quantile(data, 0.25);
            Q3 = quantile(data, 0.75);
            IQR = Q3 - Q1;
            
            lower_bound = Q1 - threshold * IQR;
            upper_bound = Q3 + threshold * IQR;
            
            outlier_flags = data < lower_bound | data > upper_bound;
            
        case 'mad'
            % Median Absolute Deviation method
            median_data = median(data);
            mad_data = median(abs(data - median_data));
            
            normalized_deviations = abs(data - median_data) / (1.4826 * mad_data);
            outlier_flags = normalized_deviations > threshold;
            
        case 'zscore'
            % Z-score method (assumes normal distribution)
            z_scores = abs(data - mean(data)) / std(data);
            outlier_flags = z_scores > threshold;
            
        case 'velocity'
            % Velocity-based outlier detection
            if length(data) < 3
                % Not enough points for velocity calculation
                data_clean = data;
                return;
            end
            
            % Compute point-to-point velocities
            dt = diff(time);
            velocities = diff(data) ./ dt;
            
            % Detect outliers in velocity space
            mad_vel = median(abs(velocities - median(velocities)));
            normalized_vel = abs(velocities - median(velocities)) / (1.4826 * mad_vel);
            
            outlier_velocities = normalized_vel > threshold;
            
            % Mark corresponding data points as outliers
            outlier_flags(2:end) = outlier_flags(2:end) | outlier_velocities;
            outlier_flags(1:end-1) = outlier_flags(1:end-1) | outlier_velocities;
    end
    
    % Create cleaned data
    data_clean = data;
    data_clean(outlier_flags) = NaN;
    
    if sum(outlier_flags) > 0
        fprintf('Detected %d outliers (%.1f%% of data)\n', ...
                sum(outlier_flags), 100*sum(outlier_flags)/n_obs);
    end
end
```

**Physical Reasonableness Checks**:
```matlab
function [data_clean, flags] = physicalQualityCheck(data, component)
    % Check for physically reasonable values
    
    flags = false(size(data));
    
    switch component
        case 'vertical'
            % Vertical displacements should be within reasonable bounds
            % Allow for ±1m total range (very conservative)
            range_check = abs(data) > 1.0; % 1 meter
            
            % Check for unrealistic velocities (>1 cm/day)
            if length(data) > 1
                velocities = abs(diff(data));
                vel_check = [false; velocities > 0.01]; % 1 cm/day
            else
                vel_check = false(size(data));
            end
            
            flags = range_check | vel_check;
            
        case 'horizontal'
            % Horizontal displacements (more variable due to tectonics)
            range_check = abs(data) > 2.0; % 2 meters
            
            if length(data) > 1
                velocities = abs(diff(data));
                vel_check = [false; velocities > 0.05]; % 5 cm/day
            else
                vel_check = false(size(data));
            end
            
            flags = range_check | vel_check;
    end
    
    data_clean = data;
    data_clean(flags) = NaN;
    
    if sum(flags) > 0
        fprintf('Failed physical checks: %d observations\n', sum(flags));
    end
end
```

### 5.4 Step 3: Gap Analysis and Data Continuity

**Objective**: Analyze temporal gaps and assess data continuity for reliable trend estimation.

**Gap Detection**:
```matlab
function gap_info = analyzeDataGaps(time, data)
    % Analyze temporal gaps in GPS time series
    
    % Remove NaN data points
    valid_idx = ~isnan(data);
    time_valid = time(valid_idx);
    
    if length(time_valid) < 2
        gap_info = struct('n_gaps', 0, 'total_gap_days', 0, 'max_gap_days', 0);
        return;
    end
    
    % Compute time differences
    dt = diff(time_valid);
    
    % Expected sampling interval (mode of time differences)
    expected_dt = mode(round(dt * 10) / 10); % Round to 0.1 days
    
    % Define gap threshold (e.g., 3x expected interval)
    gap_threshold = 3 * expected_dt;
    
    % Find gaps
    gap_idx = dt > gap_threshold;
    gap_lengths = dt(gap_idx);
    
    % Compile gap information
    gap_info = struct();
    gap_info.n_gaps = sum(gap_idx);
    gap_info.total_gap_days = sum(gap_lengths);
    gap_info.max_gap_days = max([0; gap_lengths]);
    gap_info.gap_starts = time_valid(find(gap_idx));
    gap_info.gap_ends = time_valid(find(gap_idx) + 1);
    gap_info.gap_lengths_days = gap_lengths;
    
    % Data completeness assessment
    total_span = time_valid(end) - time_valid(1);
    data_span = total_span - gap_info.total_gap_days;
    gap_info.completeness_percent = 100 * data_span / total_span;
    
    % Summary
    fprintf('Gap Analysis Summary:\n');
    fprintf('  Number of gaps: %d\n', gap_info.n_gaps);
    fprintf('  Total gap time: %.1f days\n', gap_info.total_gap_days);
    fprintf('  Maximum gap: %.1f days\n', gap_info.max_gap_days);
    fprintf('  Data completeness: %.1f%%\n', gap_info.completeness_percent);
    
    % Quality assessment
    if gap_info.completeness_percent < 70
        fprintf('  WARNING: Low data completeness (<70%%)\n');
    elseif gap_info.max_gap_days > 365
        fprintf('  WARNING: Long data gaps (>1 year) detected\n');
    else
        fprintf('  Data continuity: ACCEPTABLE\n');
    end
end
```

### 5.5 Step 4: Polynomial Detrending

**Objective**: Remove secular trends to isolate loading-related signals.

**Adaptive Detrending**:
```matlab
function [detrended_series, trend_info] = adaptiveDetrend(time, data, max_degree)
    % Adaptive polynomial detrending with optimal degree selection
    
    if nargin < 3, max_degree = 2; end
    
    % Remove NaN values for trend fitting
    valid_idx = ~isnan(data) & ~isnan(time);
    time_valid = time(valid_idx);
    data_valid = data(valid_idx);
    
    if length(data_valid) < max_degree + 2
        error('Insufficient data points for polynomial fitting');
    end
    
    % Test different polynomial degrees
    degrees = 0:max_degree;
    aic_values = zeros(size(degrees));
    bic_values = zeros(size(degrees));
    
    fprintf('Testing polynomial degrees 0-%d...\n', max_degree);
    
    for i = 1:length(degrees)
        p = degrees(i);
        
        % Fit polynomial
        [~, trend_params] = polynomialDetrend(time_valid, data_valid, p);
        
        % Compute residuals
        A = zeros(length(time_valid), p+1);
        for j = 0:p
            A(:, j+1) = (time_valid - time_valid(1)).^j;
        end
        fitted = A * trend_params;
        residuals = data_valid - fitted;
        
        % Compute information criteria
        n = length(residuals);
        k = p + 1; % Number of parameters
        rss = sum(residuals.^2); % Residual sum of squares
        
        % Akaike Information Criterion
        aic_values(i) = 2*k + n*log(rss/n);
        
        % Bayesian Information Criterion  
        bic_values(i) = k*log(n) + n*log(rss/n);
        
        fprintf('  Degree %d: AIC=%.2f, BIC=%.2f, RSS=%.6f\n', p, aic_values(i), bic_values(i), rss);
    end
    
    % Select optimal degree (minimum BIC)
    [~, optimal_idx] = min(bic_values);
    optimal_degree = degrees(optimal_idx);
    
    fprintf('Selected polynomial degree: %d\n', optimal_degree);
    
    % Fit final model with optimal degree
    [detrended_data, trend_params] = polynomialDetrend(time_valid, data_valid, optimal_degree);
    
    % Reconstruct full time series (including NaN positions)
    detrended_series = NaN(size(data));
    detrended_series(valid_idx) = detrended_data;
    
    % Create trend information structure
    trend_info = struct();
    trend_info.degree = optimal_degree;
    trend_info.parameters = trend_params;
    trend_info.aic = aic_values(optimal_idx);
    trend_info.bic = bic_values(optimal_idx);
    
    % Interpret trend parameters
    if optimal_degree >= 1
        trend_info.velocity_m_per_year = trend_params(2);
        trend_info.velocity_mm_per_year = trend_params(2) * 1000;
    end
    
    if optimal_degree >= 2
        trend_info.acceleration_m_per_year2 = 2 * trend_params(3);
        trend_info.acceleration_mm_per_year2 = 2 * trend_params(3) * 1000;
    end
    
    fprintf('Trend analysis results:\n');
    if optimal_degree >= 1
        fprintf('  Velocity: %.2f mm/year\n', trend_info.velocity_mm_per_year);
    end
    if optimal_degree >= 2
        fprintf('  Acceleration: %.2f mm/year²\n', trend_info.acceleration_mm_per_year2);
    end
end
```

### 5.6 Step 5: Harmonic Analysis (Optional)

**Objective**: Extract periodic signals for detailed seasonal analysis.

**Seasonal Component Extraction**:
```matlab
function [seasonal_info] = extractSeasonalSignals(time, data)
    % Extract annual and semi-annual components
    
    % Remove NaN values
    valid_idx = ~isnan(data) & ~isnan(time);
    time_valid = time(valid_idx);
    data_valid = data(valid_idx);
    
    if length(data_valid) < 10
        error('Insufficient data for harmonic analysis');
    end
    
    % Define frequencies of interest
    frequencies = [1.0, 2.0]; % Annual and semi-annual [cycles/year]
    freq_names = {'Annual', 'Semi-annual'};
    
    % Fit harmonic model
    [harmonic_params, harmonic_signal_valid] = fitHarmonics(time_valid, data_valid, frequencies);
    
    % Reconstruct full time series
    harmonic_signal = NaN(size(data));
    harmonic_signal(valid_idx) = harmonic_signal_valid;
    
    % Extract components
    seasonal_info = struct();
    seasonal_info.total_signal = harmonic_signal;
    seasonal_info.constant = harmonic_params.constant;
    
    for i = 1:length(frequencies)
        comp_name = lower(freq_names{i});
        
        seasonal_info.(comp_name).frequency = frequencies(i);
        seasonal_info.(comp_name).amplitude_m = harmonic_params.amplitudes(i);
        seasonal_info.(comp_name).amplitude_mm = harmonic_params.amplitudes(i) * 1000;
        seasonal_info.(comp_name).phase_rad = harmonic_params.phases(i);
        seasonal_info.(comp_name).phase_deg = harmonic_params.phases(i) * 180/pi;
        
        % Convert phase to day of year for annual component
        if frequencies(i) == 1.0
            seasonal_info.(comp_name).phase_doy = harmonic_params.phases(i) * 365.25 / (2*pi);
            if seasonal_info.(comp_name).phase_doy < 0
                seasonal_info.(comp_name).phase_doy = seasonal_info.(comp_name).phase_doy + 365.25;
            end
        end
        
        fprintf('%s component:\n', freq_names{i});
        fprintf('  Amplitude: %.2f mm\n', seasonal_info.(comp_name).amplitude_mm);
        fprintf('  Phase: %.1f° (%.1f days)\n', seasonal_info.(comp_name).phase_deg, ...
               seasonal_info.(comp_name).phase_doy);
    end
    
    % Compute residuals (detrended and de-seasoned)
    seasonal_info.residuals = data - harmonic_signal;
    residual_rms = sqrt(mean(seasonal_info.residuals(valid_idx).^2)) * 1000;
    
    fprintf('Residual RMS after harmonic removal: %.2f mm\n', residual_rms);
    
    % Variance explained by harmonic components
    total_var = var(data_valid);
    residual_var = var(seasonal_info.residuals(valid_idx));
    variance_explained = 100 * (1 - residual_var / total_var);
    
    seasonal_info.variance_explained_percent = variance_explained;
    fprintf('Variance explained by harmonic components: %.1f%%\n', variance_explained);
end
```

### 5.7 Step 6: Reference Frame Corrections

**Objective**: Apply corrections for reference frame effects and common mode errors.

**Common Mode Error Removal**:
```matlab
function [corrected_data] = removeCommonModeError(gps_database, component)
    % Remove common mode errors using network approach
    
    n_stations = gps_database.station_info.n_stations;
    corrected_data = cell(n_stations, 1);
    
    % Extract all time series for the specified component
    all_data = [];
    all_times = [];
    station_indices = [];
    
    fprintf('Computing common mode error for %s component...\n', component);
    
    for i = 1:n_stations
        ts = gps_database.time_series{i};
        if isfield(ts, component)
            data = ts.(component);
            time = ts.mjd;
            
            % Remove trend from each station first
            valid_idx = ~isnan(data);
            if sum(valid_idx) > 20 % Minimum data requirement
                [detrended, ~] = polynomialDetrend(time(valid_idx), data(valid_idx), 1);
                
                % Store for common mode computation
                all_data = [all_data; detrended];
                all_times = [all_times; time(valid_idx)];
                station_indices = [station_indices; i*ones(sum(valid_idx), 1)];
            end
        end
    end
    
    if isempty(all_data)
        error('No valid data found for common mode estimation');
    end
    
    % Compute common mode on regular time grid
    time_min = min(all_times);
    time_max = max(all_times);
    time_grid = time_min:30:time_max; % Monthly grid
    
    common_mode_signal = zeros(size(time_grid));
    
    for t = 1:length(time_grid)
        % Find data within ±15 days of this time
        time_window = abs(all_times - time_grid(t)) <= 15;
        
        if sum(time_window) >= 3 % Minimum stations requirement
            window_data = all_data(time_window);
            
            % Compute robust average (median)
            common_mode_signal(t) = median(window_data);
        else
            common_mode_signal(t) = NaN;
        end
    end
    
    % Interpolate common mode to original time stamps and subtract
    for i = 1:n_stations
        ts = gps_database.time_series{i};
        if isfield(ts, component)
            time = ts.mjd;
            data = ts.(component);
            
            % Interpolate common mode
            common_mode_interp = interp1(time_grid, common_mode_signal, time, 'linear', NaN);
            
            % Remove common mode
            corrected_data{i} = data - common_mode_interp;
        else
            corrected_data{i} = [];
        end
    end
    
    % Statistics
    cm_rms = sqrt(mean(common_mode_signal(~isnan(common_mode_signal)).^2)) * 1000;
    fprintf('Common mode RMS: %.2f mm\n', cm_rms);
end
```

### 5.8 Step 7: Uncertainty Propagation

**Objective**: Properly propagate uncertainties through the processing chain.

**Error Propagation**:
```matlab
function [uncertainties] = propagateUncertainties(gps_data, processing_info)
    % Propagate uncertainties through GPS processing steps
    
    % Extract formal uncertainties from TENV3 data
    formal_errors = gps_data.sig_up; % Standard deviations [m]
    
    % Initialize uncertainty structure
    uncertainties = struct();
    uncertainties.formal_mm = formal_errors * 1000;
    
    % Detrending uncertainty contribution
    n_trend_params = processing_info.trend_degree + 1;
    n_obs = length(gps_data.up);
    
    % Degrees of freedom loss due to trend fitting
    dof_loss_trend = n_trend_params;
    
    % Residual variance after detrending
    residuals = processing_info.detrended_residuals;
    valid_residuals = residuals(~isnan(residuals));
    
    if length(valid_residuals) > dof_loss_trend
        residual_variance = var(valid_residuals);
        trend_uncertainty_mm = sqrt(residual_variance * n_trend_params / n_obs) * 1000;
    else
        trend_uncertainty_mm = 0;
    end
    
    uncertainties.trend_contribution_mm = trend_uncertainty_mm;
    
    % Harmonic fitting uncertainty (if applicable)
    if isfield(processing_info, 'harmonic_residuals')
        n_harmonic_params = 2 * length(processing_info.frequencies) + 1;
        dof_loss_harmonic = n_harmonic_params;
        
        harmonic_residuals = processing_info.harmonic_residuals;
        valid_harmonic = harmonic_residuals(~isnan(harmonic_residuals));
        
        if length(valid_harmonic) > dof_loss_harmonic
            harmonic_variance = var(valid_harmonic);
            harmonic_uncertainty_mm = sqrt(harmonic_variance * n_harmonic_params / n_obs) * 1000;
        else
            harmonic_uncertainty_mm = 0;
        end
        
        uncertainties.harmonic_contribution_mm = harmonic_uncertainty_mm;
    else
        uncertainties.harmonic_contribution_mm = 0;
    end
    
    % Total combined uncertainty
    uncertainties.total_mm = sqrt(mean(uncertainties.formal_mm.^2) + ...
                                 trend_uncertainty_mm^2 + ...
                                 uncertainties.harmonic_contribution_mm^2);
    
    fprintf('Uncertainty budget:\n');
    fprintf('  Formal errors (RMS): %.2f mm\n', sqrt(mean(uncertainties.formal_mm.^2)));
    fprintf('  Trend fitting: ±%.2f mm\n', trend_uncertainty_mm);
    fprintf('  Harmonic fitting: ±%.2f mm\n', uncertainties.harmonic_contribution_mm);
    fprintf('  Total combined: ±%.2f mm\n', uncertainties.total_mm);
end
```

---

## 6. Statistical Methods

### 6.1 Robust Statistics

**Median Absolute Deviation (MAD)**:
```matlab
function mad_value = robustMAD(data)
    % Compute robust measure of variability
    valid_data = data(~isnan(data));
    median_data = median(valid_data);
    mad_value = median(abs(valid_data - median_data));
end
```

**Robust Correlation**:
```matlab
function [r_robust, p_value] = robustCorrelation(x, y, method)
    % Compute robust correlation coefficient
    if nargin < 3, method = 'spearman'; end
    
    % Remove pairs with missing data
    valid_pairs = ~isnan(x) & ~isnan(y);
    x_clean = x(valid_pairs);
    y_clean = y(valid_pairs);
    
    if length(x_clean) < 3
        r_robust = NaN;
        p_value = NaN;
        return;
    end
    
    switch method
        case 'spearman'
            [r_robust, p_value] = corr(x_clean, y_clean, 'type', 'Spearman');
            
        case 'kendall'
            [r_robust, p_value] = corr(x_clean, y_clean, 'type', 'Kendall');
            
        case 'robust'
            % Biweight midcorrelation (more robust than Pearson)
            r_robust = biweightMidcorr(x_clean, y_clean);
            % Approximate p-value using Fisher transformation
            n = length(x_clean);
            if n > 3
                t_stat = r_robust * sqrt((n-2)/(1-r_robust^2));
                p_value = 2 * (1 - tcdf(abs(t_stat), n-2));
            else
                p_value = NaN;
            end
    end
end
```

### 6.2 Time Series Bootstrap

**Block Bootstrap for Time Series**:
```matlab
function [bootstrap_stats] = blockBootstrap(data, n_bootstrap, block_length, statistic_func)
    % Block bootstrap for time series data to account for temporal correlation
    
    if nargin < 4
        statistic_func = @mean; % Default statistic
    end
    
    valid_data = data(~isnan(data));
    n_data = length(valid_data);
    
    if n_data < block_length
        error('Data length must be greater than block length');
    end
    
    % Pre-allocate bootstrap statistics
    bootstrap_stats = zeros(n_bootstrap, 1);
    
    for b = 1:n_bootstrap
        % Generate bootstrap sample using overlapping blocks
        n_blocks = ceil(n_data / block_length);
        bootstrap_sample = zeros(n_blocks * block_length, 1);
        
        for i = 1:n_blocks
            % Random starting position for block
            start_pos = randi(n_data - block_length + 1);
            block_data = valid_data(start_pos:start_pos + block_length - 1);
            
            % Insert block into bootstrap sample
            insert_start = (i-1) * block_length + 1;
            insert_end = insert_start + block_length - 1;
            bootstrap_sample(insert_start:insert_end) = block_data;
        end
        
        % Trim to original data length
        bootstrap_sample = bootstrap_sample(1:n_data);
        
        % Compute statistic for this bootstrap sample
        bootstrap_stats(b) = statistic_func(bootstrap_sample);
    end
end
```

---

## 7. Quality Control and Validation

### 7.1 Data Quality Metrics

**Completeness Assessment**:
```matlab
function [quality_metrics] = assessDataQuality(gps_data)
    % Comprehensive data quality assessment
    
    quality_metrics = struct();
    
    % Basic statistics
    total_obs = length(gps_data.up);
    valid_obs = sum(~isnan(gps_data.up));
    quality_metrics.completeness_percent = 100 * valid_obs / total_obs;
    
    % Temporal coverage
    time_span_years = (max(gps_data.mjd) - min(gps_data.mjd)) / 365.25;
    quality_metrics.time_span_years = time_span_years;
    
    % Sampling rate analysis
    valid_times = gps_data.mjd(~isnan(gps_data.up));
    if length(valid_times) > 1
        time_diffs = diff(valid_times);
        quality_metrics.median_sampling_days = median(time_diffs);
        quality_metrics.sampling_regularity = std(time_diffs) / mean(time_diffs);
    end
    
    % Precision indicators
    valid_data = gps_data.up(~isnan(gps_data.up));
    if ~isempty(valid_data)
        quality_metrics.data_scatter_mm = std(valid_data) * 1000;
        quality_metrics.data_range_mm = (max(valid_data) - min(valid_data)) * 1000;
    end
    
    % Formal uncertainty statistics
    if isfield(gps_data, 'sig_up')
        valid_errors = gps_data.sig_up(~isnan(gps_data.sig_up));
        if ~isempty(valid_errors)
            quality_metrics.mean_formal_error_mm = mean(valid_errors) * 1000;
            quality_metrics.formal_error_variability = std(valid_errors) / mean(valid_errors);
        end
    end
    
    % Quality classification
    if quality_metrics.completeness_percent >= 80 && time_span_years >= 5
        quality_metrics.overall_grade = 'EXCELLENT';
    elseif quality_metrics.completeness_percent >= 60 && time_span_years >= 3
        quality_metrics.overall_grade = 'GOOD';
    elseif quality_metrics.completeness_percent >= 40 && time_span_years >= 2
        quality_metrics.overall_grade = 'FAIR';
    else
        quality_metrics.overall_grade = 'POOR';
    end
    
    % Display summary
    fprintf('Data Quality Assessment:\n');
    fprintf('  Completeness: %.1f%% (%d/%d observations)\n', ...
            quality_metrics.completeness_percent, valid_obs, total_obs);
    fprintf('  Time span: %.1f years\n', time_span_years);
    fprintf('  Data scatter: %.2f mm\n', quality_metrics.data_scatter_mm);
    fprintf('  Overall grade: %s\n', quality_metrics.overall_grade);
end
```

### 7.2 Physical Validation

**Expected Signal Characteristics**:
```matlab
function [validation_results] = validatePhysicalBehavior(processed_data, station_info)
    % Validate GPS signals against expected physical behavior
    
    validation_results = struct();
    validation_results.tests_passed = 0;
    validation_results.total_tests = 5;
    
    % Test 1: Seasonal amplitude range
    if isfield(processed_data, 'seasonal_info')
        annual_amplitude_mm = processed_data.seasonal_info.annual.amplitude_mm;
        
        % Expected range: 1-20 mm for hydrological loading
        if annual_amplitude_mm >= 1 && annual_amplitude_mm <= 20
            validation_results.amplitude_test = 'PASS';
            validation_results.tests_passed = validation_results.tests_passed + 1;
        else
            validation_results.amplitude_test = 'FAIL';
        end
        
        fprintf('Seasonal amplitude test: %.2f mm - %s\n', ...
                annual_amplitude_mm, validation_results.amplitude_test);
    else
        validation_results.amplitude_test = 'SKIPPED';
    end
    
    % Test 2: Phase consistency
    if isfield(processed_data, 'seasonal_info')
        phase_doy = processed_data.seasonal_info.annual.phase_doy;
        
        % Expected: maximum loading in winter/spring (day 30-150)
        % GPS shows maximum UPLIFT when loading decreases
        if (phase_doy >= 150 && phase_doy <= 300) % Summer/fall uplift
            validation_results.phase_test = 'PASS';
            validation_results.tests_passed = validation_results.tests_passed + 1;
        else
            validation_results.phase_test = 'QUESTIONABLE';
        end
        
        fprintf('Phase timing test: Day %.0f - %s\n', phase_doy, validation_results.phase_test);
    else
        validation_results.phase_test = 'SKIPPED';
    end
    
    % Test 3: Velocity reasonableness
    if isfield(processed_data, 'trend_info') && isfield(processed_data.trend_info, 'velocity_mm_per_year')
        velocity_mm_yr = processed_data.trend_info.velocity_mm_per_year;
        
        % Most GPS sites have secular velocities < 50 mm/year
        if abs(velocity_mm_yr) <= 50
            validation_results.velocity_test = 'PASS';
            validation_results.tests_passed = validation_results.tests_passed + 1;
        else
            validation_results.velocity_test = 'QUESTIONABLE';
        end
        
        fprintf('Secular velocity test: %.2f mm/yr - %s\n', velocity_mm_yr, validation_results.velocity_test);
    else
        validation_results.velocity_test = 'SKIPPED';
    end
    
    % Test 4: Noise level
    if isfield(processed_data, 'residuals')
        residual_rms_mm = sqrt(mean(processed_data.residuals(~isnan(processed_data.residuals)).^2)) * 1000;
        
        % Typical GPS noise: 1-5 mm RMS
        if residual_rms_mm >= 0.5 && residual_rms_mm <= 10
            validation_results.noise_test = 'PASS';
            validation_results.tests_passed = validation_results.tests_passed + 1;
        else
            validation_results.noise_test = 'QUESTIONABLE';
        end
        
        fprintf('Noise level test: %.2f mm RMS - %s\n', residual_rms_mm, validation_results.noise_test);
    else
        validation_results.noise_test = 'SKIPPED';
    end
    
    % Test 5: Geographic consistency
    latitude = station_info.latitude;
    
    if isfield(processed_data, 'seasonal_info')
        amplitude_mm = processed_data.seasonal_info.annual.amplitude_mm;
        
        % Higher latitudes often show larger hydrological signals
        if latitude > 40 % High latitude
            expected_range = [2, 15]; % mm
        else % Lower latitude
            expected_range = [1, 8];  % mm
        end
        
        if amplitude_mm >= expected_range(1) && amplitude_mm <= expected_range(2)
            validation_results.geographic_test = 'PASS';
            validation_results.tests_passed = validation_results.tests_passed + 1;
        else
            validation_results.geographic_test = 'QUESTIONABLE';
        end
        
        fprintf('Geographic consistency: %.2f mm at %.1f°N - %s\n', ...
                amplitude_mm, latitude, validation_results.geographic_test);
    else
        validation_results.geographic_test = 'SKIPPED';
    end
    
    % Overall validation result
    validation_results.pass_rate = validation_results.tests_passed / validation_results.total_tests;
    
    if validation_results.pass_rate >= 0.8
        validation_results.overall_validation = 'EXCELLENT';
    elseif validation_results.pass_rate >= 0.6
        validation_results.overall_validation = 'GOOD';
    elseif validation_results.pass_rate >= 0.4
        validation_results.overall_validation = 'ACCEPTABLE';
    else
        validation_results.overall_validation = 'POOR';
    end
    
    fprintf('\nValidation Summary: %d/%d tests passed (%.0f%%) - %s\n', ...
            validation_results.tests_passed, validation_results.total_tests, ...
            100*validation_results.pass_rate, validation_results.overall_validation);
end
```

---

## 8. Error Analysis

### 8.1 Error Budget

**GPS Error Sources**:

1. **Observation Errors** (~1-3 mm):
   - Receiver noise
   - Multipath effects
   - Atmospheric delays

2. **Processing Errors** (~0.5-2 mm):
   - Reference frame realization
   - Antenna calibration
   - Tide model errors

3. **Analysis Errors** (~0.5-1 mm):
   - Trend model selection
   - Outlier detection sensitivity
   - Temporal interpolation

**Error Combination**:
```matlab
function [total_uncertainty] = combineErrorSources(error_components)
    % Combine independent error sources using root-sum-squares
    
    % Error components in mm
    observation_error = error_components.observation_mm;
    processing_error = error_components.processing_mm;
    analysis_error = error_components.analysis_mm;
    
    % Total uncertainty (assuming independence)
    total_uncertainty = sqrt(observation_error^2 + processing_error^2 + analysis_error^2);
    
    fprintf('GPS Error Budget:\n');
    fprintf('  Observation errors: ±%.2f mm\n', observation_error);
    fprintf('  Processing errors:  ±%.2f mm\n', processing_error);
    fprintf('  Analysis errors:    ±%.2f mm\n', analysis_error);
    fprintf('  Total uncertainty:  ±%.2f mm\n', total_uncertainty);
    
    % Individual contributions
    obs_contribution = 100 * observation_error^2 / total_uncertainty^2;
    proc_contribution = 100 * processing_error^2 / total_uncertainty^2;
    analysis_contribution = 100 * analysis_error^2 / total_uncertainty^2;
    
    fprintf('Error contributions:\n');
    fprintf('  Observation: %.1f%%\n', obs_contribution);
    fprintf('  Processing:  %.1f%%\n', proc_contribution);
    fprintf('  Analysis:    %.1f%%\n', analysis_contribution);
end
```

### 8.2 Temporal Correlation Effects

**Allan Variance Analysis**:
```matlab
function [tau, allan_var] = computeAllanVariance(data, sampling_interval)
    % Compute Allan variance to characterize temporal correlations
    
    % Remove NaN values and ensure uniform sampling
    valid_data = data(~isnan(data));
    n = length(valid_data);
    
    % Time scales to evaluate
    max_tau = floor(n/3); % Maximum averaging time
    tau_indices = unique(round(logspace(0, log10(max_tau), 20)));
    
    allan_var = zeros(size(tau_indices));
    tau = tau_indices * sampling_interval;
    
    for i = 1:length(tau_indices)
        m = tau_indices(i);
        
        if m >= n/2
            allan_var(i) = NaN;
            continue;
        end
        
        % Compute overlapping Allan variance
        sum_var = 0;
        count = 0;
        
        for j = 1:(n - 2*m)
            y1 = mean(valid_data(j:j+m-1));
            y2 = mean(valid_data(j+m:j+2*m-1));
            sum_var = sum_var + (y2 - y1)^2;
            count = count + 1;
        end
        
        if count > 0
            allan_var(i) = sum_var / (2 * count);
        else
            allan_var(i) = NaN;
        end
    end
    
    % Remove NaN values
    valid_idx = ~isnan(allan_var);
    tau = tau(valid_idx);
    allan_var = allan_var(valid_idx);
end
```

---

## Summary

This Part 2 documentation provides comprehensive coverage of GPS time series analysis for crustal deformation studies. The methodology enables extraction of hydrological loading signals from GPS vertical displacement measurements with typical uncertainties of ±1-3 mm.

**Key Accomplishments**:

1. **Complete Processing Pipeline**: From raw TENV3 files to analysis-ready time series
2. **Robust Statistical Methods**: Outlier detection, gap analysis, and uncertainty quantification  
3. **Physical Validation**: Tests against expected loading signal characteristics
4. **Quality Control Framework**: Comprehensive data quality assessment and validation

The processed GPS time series provide the observational foundation for comparison with GRACE-based theoretical predictions in Part 3, enabling validation of elastic loading theory and assessment of GRACE data quality.

**Typical Processing Results**:
- Secular velocities: ±0.1-50 mm/year (tectonic + GIA)
- Annual amplitudes: 1-20 mm (hydrological loading)
- Residual noise: 1-5 mm RMS (instrument + processing)
- Total uncertainty: ±2-5 mm (combined error sources)

---

*This documentation provides complete theoretical background and practical implementation guidance for GPS time series analysis without requiring external references beyond the foundational geodetic literature.*