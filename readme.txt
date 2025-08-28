Title: Modeling the elastic response of Earth’s crust to hydrological loading using GRACE data and
comparison with GPS observations
Summary: The Earth deforms elastically under gravitational force of surface loads including
atmospheric, oceanic and hydrological masses. The hydrological loading is often the largest component
of surface mass changes. GRACE satellites provide monthly variations of total water storage at scales
of 300 – 500 km. Using the elastic loading theory, these forcing data can be used to model crustal
deformation. As an independent technique, GPS is capable of measuring crustal deformation resulting –
among many other processes – from hydrological loading. In this project, the student will use monthly
GRACE spherical harmonic coefficients from 2002 to 2017 and model the vertical deformation using a
PREM Earth model. The modeled deformation are compared with GPS height time series at GPS sites
in a study area.


Data for this project is available here in Data/

Spherical harmonic functions are in Satellite_geodesy_assignment/
More helping functino at /project/functinos

template.m to mimic and then filled to finsih the assignement. 

  Part 1 - GRACE to Vertical Deformation (Spherical Harmonic Synthesis):
  - Load GRACE spherical harmonic coefficients
  - Apply corrections (C20 replacement, degree-1 coefficients)
  - Load Love numbers for PREM model
  - Set up spatial grid (colatitude/longitude)
  - Compute associated Legendre functions
  - Apply spherical harmonic synthesis formula
  - Convert to physical deformation using scaling factors
  - Output: vertical deformation field

  Part 2 - GPS Height Time Series Analysis:
  - Load GPS station coordinates and time series
  - Remove secular trends (polynomial fitting)
  - Remove seasonal signals if needed
  - Handle outliers and data gaps
  - Extract vertical displacement variations
  - Convert to same reference frame as GRACE
  - Output: GPS vertical displacements

  Part 3 - Comparison and Validation:
  - Spatial interpolation of GRACE results to GPS locations
  - Temporal alignment of both datasets
  - Statistical comparison (correlation, RMSE, etc.)
  - Seasonal analysis and phase comparisons
  - Identify systematic differences
  - Validate theoretical predictions