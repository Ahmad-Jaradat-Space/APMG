function constants = physicalConstants()
% physicalConstants - Centralized physical constants for GRACE analysis
%
% SYNTAX:
%   constants = physicalConstants()
%
% OUTPUT:
%   constants - Structure containing physical constants used in geodesy
%
% CONSTANTS:
%   GM           - Gravitational parameter [m³/s²]
%   R            - Earth mean radius [m]
%   rho_water    - Water density [kg/m³]
%   rho_earth    - Mean Earth density [kg/m³]
%   g            - Mean surface gravity [m/s²]
%   omega        - Earth rotation rate [rad/s]
%   flattening   - Earth flattening factor
%   a            - Semi-major axis [m]
%   b            - Semi-minor axis [m]
%
% REFERENCES:
%   - IERS Conventions 2010
%   - GRS80 Reference Ellipsoid
%   - Consistent with main_fun_new.m implementation
%
% Author: GRACE Analysis Project
% Date: 2025

% Initialize constants structure
constants = struct();

%% Fundamental Constants (IERS Conventions 2010)
constants.GM = 3.986004418e14;          % Gravitational parameter [m³/s²]
constants.R = 6378137.0;                % Earth mean radius [m] (GRS80)
constants.omega = 7.2921159e-5;         % Earth rotation rate [rad/s]

%% Density Constants
constants.rho_water = 1000.0;           % Water density [kg/m³]
constants.rho_earth = 5517.0;           % Mean Earth density [kg/m³]
constants.rho_ice = 917.0;              % Ice density [kg/m³]

%% Reference Ellipsoid (GRS80)
constants.a = 6378137.0;                % Semi-major axis [m]
constants.f = 1/298.257222101;          % Flattening
constants.b = constants.a * (1 - constants.f);  % Semi-minor axis [m]

%% Derived Constants
constants.g = 9.80665;                  % Standard gravity [m/s²]
constants.c = 299792458;                % Speed of light [m/s]

%% GRACE-specific Constants
constants.grace_altitude = 450000;      % Typical GRACE altitude [m]
constants.grace_inclination = 89;       % GRACE inclination [degrees]

%% Unit Conversions
constants.deg2rad = pi/180;             % Degrees to radians
constants.rad2deg = 180/pi;             % Radians to degrees
constants.m2mm = 1000;                  % Meters to millimeters
constants.mm2m = 1/1000;                % Millimeters to meters
constants.mgal2ms2 = 1e-5;              % Milligal to m/s²
constants.ms2_2mgal = 1e5;              % m/s² to milligal

%% Time Constants
constants.mjd_j2000 = 51544.5;          % MJD of J2000.0 epoch
constants.days_per_year = 365.25;       % Days per year (Julian)
constants.seconds_per_day = 86400;      % Seconds per day

%% Love Numbers Reference Values (PREM)
% Load Love numbers from Wang et al. (2012) - CORRECTED POSITIVE VALUES
constants.h2_prem = 0.6149;            % Vertical displacement Love number
constants.l2_prem = 0.0839;            % Horizontal displacement Love number  
constants.k2_prem = 0.3020;            % Gravitational potential Love number

%% Scaling Factors for Different Quantities
% Based on Wahr et al. (1998) and main_fun_new.m implementation

% Disturbance potential factor
constants.dist_potential_factor = constants.GM / constants.R;

% Geoid height factor  
constants.geoid_factor = constants.R;

% Gravity anomaly factor (function of degree n)
% anomaly_factor(n) = GM/R² * (n+1)
constants.gravity_anomaly_base = constants.GM / constants.R^2;

% Height factor for satellite altitude (function of degree n and height h)
% height_factor(n,h) = (R/(R+h))^n
constants.height_factor_base = @(h) constants.R / (constants.R + h);

% Vertical deformation scaling (Wahr et al. 1998) 
% Exact formula from Wahr et al. (1998): (R × ρ_water) / (3 × ρ_earth)
% This converts dimensionless spherical harmonic coefficients to meters of displacement
constants.wahr_scaling_factor = (constants.R * constants.rho_water) / (3 * constants.rho_earth);  % [m]

%% Quality Control
% Validate critical constants
if constants.GM <= 0
    error('Invalid GM value');
end

if constants.R <= 0 || constants.R > 1e7
    error('Invalid Earth radius');
end

if constants.rho_water <= 0 || constants.rho_earth <= 0
    error('Invalid density values');
end

%% Display Summary
fprintf('Physical constants loaded:\n');
fprintf('  GM = %.6e m³/s²\n', constants.GM);
fprintf('  R = %.1f m\n', constants.R);
fprintf('  ρ_water = %.1f kg/m³\n', constants.rho_water);
fprintf('  ρ_earth = %.1f kg/m³\n', constants.rho_earth);
fprintf('  Wahr et al. (1998) scaling factor = %.1f m\n', constants.wahr_scaling_factor);

end