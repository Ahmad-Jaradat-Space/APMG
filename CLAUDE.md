# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a MATLAB-based geodetic data analysis project for modeling the elastic response of Earth's crust to hydrological loading using GRACE satellite data and GPS observations. The project is part of Advanced Physical and Mathematical Geodesy (APMG) coursework.

## Project Structure

```
APMG/
├── Satellite_geodesy_assignment/  # Satellite altimetry analysis code
│   ├── lib/matlab/                # Library functions for altimetry
│   └── src/                       # Source code
├── project/                       # Main GRACE project
│   ├── data/                      
│   │   ├── grace/                 # GRACE spherical harmonic coefficients (.gfc files)
│   │   ├── gps/                   # GPS station time series (.tenv3 files)
│   │   └── aux/                   # Auxiliary data (C20 coefficients, degree-1 coefficients)
│   ├── functions/                 # Core analysis functions
│   └── output/                    # Results directory
└── template.m                     # Main template for analysis

```

## Key Commands

### Running MATLAB Scripts
```matlab
% To run the main analysis
cd project/
matlab -r "template"

% To add all subdirectories to path
addpath(genpath(pwd));
```

### Data Processing Workflow
1. Load GPS time series data using `load_tenv3()`
2. Read GRACE spherical harmonic coefficients using `readSHC()`
3. Compute Legendre functions using `pnm()`
4. Fit polynomial trends using `fitPolynomial()`
5. Calculate elastic deformation from GRACE data
6. Compare with GPS observations

## Core Functions

### Data I/O Functions
- `load_tenv3(filename)` - Loads GPS time series in tenv3 format
- `readSHC(file)` - Reads spherical harmonic coefficients from GRACE files
- `readConfig(configfile)` - Reads configuration files

### Mathematical Functions
- `pnm(n, theta, dm)` - Computes fully normalized associated Legendre functions
- `fitPolynomial(t, y, p, sigma0)` - Fits polynomial of degree p using least squares
- `solve_LSA()` - Solves least squares adjustment problems
- `get_designmatrx_harmonic()` - Creates design matrix for harmonic analysis

### Analysis Functions
- `outlier_test()` - Statistical outlier detection
- `significance_test()` - Statistical significance testing
- `globalTest()` - Global model test
- `select_earthquakes()` - Earthquake event selection for GPS data

### Time Conversion Functions
- Various MJD (Modified Julian Date) conversion utilities in `nsea_project/`

## Data Formats

### GPS Data (tenv3)
- Columns: time (MJD), east, north, up components with uncertainties
- Files contain station positions and displacements

### GRACE Data (.gfc)
- Spherical harmonic coefficients (Cnm, Snm)
- DDK3 filtered data from CSR (UTCSR)
- Time period: 2002-2017
- Format: degree (n), order (m), Cnm, Snm coefficients

## Analysis Parameters

### Default Settings
- Maximum spherical harmonic degree: typically 60
- Polynomial trend degree: 2
- Earth model: PREM (Preliminary Reference Earth Model)
- Coordinate system: WGS84

## Important Notes

1. **Memory Management**: GRACE data processing is memory-intensive. Clear variables when not needed.

2. **Path Management**: Always use `addpath(genpath(pwd))` to ensure all functions are accessible.

3. **Time Systems**: The project uses Modified Julian Date (MJD) for time representation. Conversion functions are available in `nsea_project/`.

4. **Filtering**: GRACE data is pre-filtered with DDK3 filter to reduce noise.

5. **Covariance Matrices**: Functions handle full covariance matrices for proper error propagation.

## Key Physical Concepts

- **Elastic Loading**: Earth's crust deforms elastically under surface mass loads
- **Spherical Harmonics**: Mathematical representation of Earth's gravity field
- **Hydrological Loading**: Dominant component of seasonal crustal deformation
- **GRACE Mission**: Measures monthly variations in Earth's gravity field

## Reference Publications

Key papers for methodology:
- Wahr et al. (1998) - Theoretical framework for loading deformation
- Fu & Freymueller (2012) - GPS-GRACE comparison methodology
- Karegar et al. (2018) - Hybrid estimation methods