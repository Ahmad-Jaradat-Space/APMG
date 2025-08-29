# Text Extraction from Provided PDFs

> Generated automatically. Page-level text was extracted directly with no paraphrasing; formatting and figure text may be incomplete due to PDF structure.

## lec-01.pdf

_Total pages: 44_

### Page 1

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Introduction to GRACE/FO and working with
their data
Makan Karegar
April 14, 2022
karegar@uni-bonn.de

### Page 2

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Gravity Recovery & Climate Experiment (GRACE)
Joint NASA and DLR mission launched in March 2002 and terminated in
October 2017 (GRACE-FO since May 2018) with a goal to measure time
variable gravity field due to changes in land water storage and ocean
circulations and movements.

### Page 3

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Concepts
High-low SST
low-low SST
SST: Satellite-to-Satellite Tracking

### Page 4

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
Concepts

### Page 5

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
Concepts

### Page 6

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Mass transport and distribution in the Earth system
Mass density is changing at different time and spatial scales.Newton’s integral for
gravitational potential:
Any changes in mass density
results in potential or gravity
Changes (eq. 1-12 in Hofmann &
Moritz, 2006)
4D density distribution
G Newtonian gravitational constant
dv an element of volume
l distance between x and x’


### Page 7

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Surface density change to potential change
Most of mass variability concentrated within a thin layer of thickness H  (~
10 km) near the surface of Earth containing atmosphere, oceans, ice
sheets, and land water storage . This layer is subject to significant mass
fluctuations.
Density of water 1000 kg/m3Let’s expand surface density over sphere with SH coefficients of         &


And expand the reciprocal distance (1/ l) to Legendre’s polynomial ( Hofmann &
Moritz 2006, Section 1.11, eq. 1-108 )
After putting SH expression of 1/ l and SH formula of surface density in the
Newton’s integral, we get the potential variation  in terms of SHCs of surface
density          &
Pn can be expanded to SH (eq. 1-105)
or Lecture Note 2, P 21


### Page 8

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
Loading effect
* Any change in surface density will cause deformation  within the solid Earth,
  leading to an additional surface density change.
* The gravity change caused by these “solid Earth mass anomalies” is up to
   30% of the gravity caused by the surface mass itself.
* It can be represented in terms of load Love numbers: kl
satellite gravimetry measures mass effect + its elastic response
Surface density change to potential change: effect of
solid Earth deformation

### Page 9

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Thin layer assumption
changes in potential in terms of SHCs of surface density
r2R: mean radius of the Earth,
H: thickness of thin leyer where most of mass changes occur (~ 10 km)
r = R and  rE = R + H
binomial expansion: [ rE/r ]n = [ ( R + H)/R ]n   =
 1 + H/R + 2*H/R + ...+ nmax*H/R

Suppose H is thin enough so [ rE/r ]n = 1
** max committing error 10/6378 = 0.001

### Page 10

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
Relationship between potential and surface density
SHCs
More in lecture 4Potential (geoid height) SHCs  ↔ Mass density SHCs
 What GRACE delivers


### Page 11

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO           Equivalent water thickness (total water storage)
Basic Earth parameters
(radius, avg. density, Love Load
numbers)Time-Variable Gravity
Coefficients
from GRACE ProjectThe dominant mass changes are related to movement of water in the oceans,
on land and through the atmosphere. Therefore one can scale by water
density (ρW) to get equivalent water thickness  (EWT, or total water storage )
change over land and ocean and compute the water storage term over large-
scales:
    EWT = surface density / pw    >>> EWT in m SI unit


### Page 12

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO GRACE will measure all gravitational variations, even those that have
time scales much shorter than one month.
In order to get the monthly variation correct, we have to model and
remove
* Ocean tides * Solid Earth tides
* Atmospheric mass changes (namely tidal and non-tidal atmospheric
loading)
* Post glacial rebound (GIA) – for contemporary mass change

GRACE estimates are non-tidal ocean mass variations, land hydrology,
and cryospheric changes.Dealiasing GRACE observations or background models
https://www.gfz-potsdam.de/en/aod1b/Atmosphere and Ocean De-aliasing Level-1B (AOD1B) products from GFZ.
AOD1B products are 3-hourly series of spherical harmonic coefficients up to
 degree and order 180 at:

### Page 13

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Commission error vs. omission error
■□□□□□
Spectrum of mean and time-variable gravity field
40 km200 km* GRACE is most sensitive to
long-wavelengths and has higher
random error at shorter
wavelengths.
* The inherent GRACE
resolution is caused by the
height of the satellites above
the Earth (~300-400 km) and
the separation between them
(~200 km) Spatial resolution and
max of SH expansion:
L ~ 20,000 km / nmax

### Page 14

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Global water cycle
© 2020 The Global Water Cycle. All Rights
Reserved.

### Page 15

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Sverdrup
Durack (2015)sverdrup (Sv): 1  million m3/s

### Page 16

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Water storage
On land, water storage changes due to:
* surface water (rivers, lakes, reservoirs flooding)
* soil moisture
* snow
* ground water
On ocean, changes due to:
* mass fluxes in and out
* Circulation changes
   due to winds
[Trenberth et al., 2007].

### Page 17

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO GRACE data: three levels of data
Level 1
•Ranging + GPS data
•Very difficult to work with
Level 2
•Spherical harmonic coefficients (SHCs) that
describe monthly gravity potential,
worldwide
•Difficult for non-specialists but doable
•Contain all information
Level 3
•Monthly maps of water storage anomalies
(1x1 deg)
•Easy to work with for non-specialists
•Does not contain all info

### Page 18

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO GRACE data - SH coefficients
Data centers provide each month (since ~4/2002) gravity SHCs(what you have to do when working with the SHCs from the data centers)

### Page 19

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO GRACE data: SH coefficients
… 3.986004415e14
  6.3781363e+06
  360
    0    0  1.0000000000000000e+00  0.0000000000000000e+00
    1    0  0.0000000000000000e+00  0.0000000000000000e+00
    1    1  0.0000000000000000e+00  0.0000000000000000e+00
    2    0 -4.8416537173600000e-04  0.0000000000000000e+00
    2    1 -1.8698763595500000e-10  1.1952801203099999e-09
    2    2  2.4391435239799998e-06 -1.4001668365399999e-06
    3    0  9.5725417379199997e-07  0.0000000000000000e+00
    3    1  2.0299888218400001e-06  2.4851315871599998e-07
    3    2  9.0462776860499998e-07 -6.1902594420499998e-07
    3    3  7.2107265705700001e-07  1.4143562695799999e-06
    4    0  5.3987386378900003e-07  0.0000000000000000e+00
    4    1 -5.3632161697100004e-07 -4.7344026585300002e-07
    4    2  3.5069410578500000e-07  6.6267157254000005e-07
    4    3  9.9077180382900004e-07 -2.0092836917700000e-07
    4    4 -1.8856080273499999e-07  3.0885316933300002e-07
    5    0  6.8532347563000006e-08  0.0000000000000000e+00
    5    1 -6.2101212852799994e-08 -9.4422612752500001e-08
    5    2  6.5243829761200005e-07 -3.2334961266800000e-07
    5    3 -4.5195540607099999e-07 -2.1484719062400001e-07
    5    4 -2.9530164765400001e-07  4.9665887676899997e-08<?xml version="1.0" encoding="UTF-8" ?>
<EarthGravityModel>
<GM>3.98600440000000e+14</GM>
<R>6.37813700000000e+06</R>
<maxDegree>360</maxDegree>
<cnm degree="0" order="0">1.00000000000000e+00</cnm>
<cnm degree="2" order="0">-4.84165532804000e-04</cnm>
<cnm degree="2" order="1">8.57179552165000e-13</cnm>
<snm degree="2" order="1">2.89607376372000e-12</snm>
<cnm degree="2" order="2">2.43815798120000e-06</cnm>
<snm degree="2" order="2">-1.39990174643000e-06</snm>
<cnm degree="3" order="0">9.57139401177000e-07</cnm>
<cnm degree="3" order="1">2.02968777310000e-06</cnm>
<snm degree="3" order="1">2.49431310090000e-07</snm>
<cnm degree="3" order="2">9.04648670700000e-07</cnm>Normalized SH coefficients

### Page 20

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Filtering
 Unfiltered GRACE data looks very noise: example here –
 total water storage


### Page 21

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Spatial Averaging to Improve Accuracy
Filtering: two steps – spatial averaging
radius of spatial filteringFor example your kernel could be:

### Page 22

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
GRACE‐derived maps of monthly anomaly of water storage, smoothed
with Gaussian (isotropic)  filters of different radius. (a) Unsmoothed;
(b) 250 km radius; (c) 500 km radius; (d) 750 km radius. (Swenson &
Wahr, 2007)Filtering: two steps – spatial average or smoothing
Isotropic Gaussian
             filter

### Page 23

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Filtering: two steps – decorrelation/Destriping filter
Monthly SHCs are correlated due
to sensor errors, orbit geometry,
and errors in background model.
Swenson & Wahr suggested to
smooth the SHCs for a particular
order (m) with a quadratic
polynomial

### Page 24

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
GRACE‐derived maps of monthly anomaly of water storage. (a) Unfiltered, no smoothing; (b)
filtered with correlated‐error filter, no smoothing; (c) unfiltered and smoothed with 500 km
Gaussian; (d) filtered with correlated‐error filter and smoothed with 500 km Gaussian.
(Swenson & Wahr, 2007)Filtering: two steps – decorrelation + smoothing

### Page 25

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
Filtering: two steps – DDK filter
Kusche designed a filter for
GRACE SHCs using the a
prior information about the
unfiltered coefficients, that is
GRACE error covariance
matrix and signal covariance
estimated from geophysical
models

### Page 26

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO GRACE early results (Tapley et al., 2004, Science)
   Orbit pattern (artefacts)
  filtering (500-800 km)

### Page 27

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO GRACE early results (Tapley et al., 2004, Science)
  EWH ~ 20 * geoid heightGRACE has revealed that global models often
underestimate the seasonal amplitude of total water storage

### Page 28

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Issues with GRACE data
•Accuracy varies over mission lifetime (due to instrument aging and
orbital pattern)
•For some months, data is missing (due to instrument problems)
•Level 2 data (SHCs): several corrections are required, some
depend on application
•GRACE data has increasing noise at shorter scales → filtering of
SHCs (level 2 data) is required
•Working with level 2 data may require seeking aid
•Level 2 data (SHCs) provided by 3 official data centres NASA JPL,
CSR U Texas, GFZ Potsdam (GRACE SDS Science Data System)
•several unofficial solutions (CNES/GRDS, Universities) available
•Several alternative data products available: masscon solutions,
daily solutions, …
•Level 3 data (total water mass grids) are easy to use

### Page 29

29
GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
Agriculture
(CA drought)Alaska GlacierNegativePositive Mass Balance
Pumping and irrigation
Earthquake
Reservoir
PrecipitationMass LossENSOAgricultureApplications: Land hydrology / water storage
direct water use +
climate change +
climate variability
Eicker et al. (2016)

### Page 30

30
GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
Trends of water mass in the GRACE era (since 2003)
Werth et al., 2018:

### Page 31

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Assimilation of GRACE data into hydrological models
•e.g. for drought monitoring
Integration of GRACE data  into models by data
assimilation

### Page 32

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
Data assimilation (DA) minimizes
Hydrological model:
•storage (soil layers,
surface water,
groundwater…)
•grid cells
•E.g. daily time steps
•Calibration parameters
Döll et al (2003):
GRACE: total water storage
anomalies (TWSA) Integration of GRACE data  into models by data
assimilation

### Page 33

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
DA minimizes
GRACEModel runs Integration of GRACE data  into models by data
assimilation
Ensemble Kalman FilterGRACE: Error variance
covariance matrix for gridded
EWH

### Page 34

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
(Eicker et al., 2014)GRACE
WGHM
DA Integration of GRACE data  into models by data
assimilation
                        2005

### Page 35

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO
 Integration of GRACE data and hydro model for Solid Earth
deformation
(Karegar et al., 2018)
Hydrological model provides TWS with
higher spatial resolutions than GRACE. A
hybrid approach was suggested to combine
GRACE and hydro model in order to better
estimate hydrological loading deformation
(more in the next lecture)

### Page 36

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Greenland ice mass imbalance due to global warming


### Page 37

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Data products 2 & 3: ICGEM
 3.986004415e14
  6.3781363e+06
  360
    0    0  1.0000000000000000e+00  0.0000000000000000e+00
    1    0  0.0000000000000000e+00  0.0000000000000000e+00
    1    1  0.0000000000000000e+00  0.0000000000000000e+00
    2    0 -4.8416537173600000e-04  0.0000000000000000e+00
    2    1 -1.8698763595500000e-10  1.1952801203099999e-09
    2    2  2.4391435239799998e-06 -1.4001668365399999e-06
    3    0  9.5725417379199997e-07  0.0000000000000000e+00
    3    1  2.0299888218400001e-06  2.4851315871599998e-07
    3    2  9.0462776860499998e-07 -6.1902594420499998e-07
    3    3  7.2107265705700001e-07  1.4143562695799999e-06
    4    0  5.3987386378900003e-07  0.0000000000000000e+00
    4    1 -5.3632161697100004e-07 -4.7344026585300002e-07
    4    2  3.5069410578500000e-07  6.6267157254000005e-07
    4    3  9.9077180382900004e-07 -2.0092836917700000e-07
    4    4 -1.8856080273499999e-07  3.0885316933300002e-07
    5    0  6.8532347563000006e-08  0.0000000000000000e+00
    5    1 -6.2101212852799994e-08 -9.4422612752500001e-08
    5    2  6.5243829761200005e-07 -3.2334961266800000e-07
    5    3 -4.5195540607099999e-07 -2.1484719062400001e-07
    5    4 -2.9530164765400001e-07  4.9665887676899997e-08Comprehensive archive:
http://icgem.gfz-potsdam.de/home

### Page 38

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO


### Page 39

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO


### Page 40

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Data products II: TWS maps
Begin-end: yyyydoy-yyyydoy
Data center /
version
Physical quantity: here TWS

### Page 41

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Data products II: TWS maps


### Page 42

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Data products 2 & 3: GRACE Tellus
•https://grace.jpl.nasa.gov/
•1*1 degree grids of TWS anomaly (with respect to long-term mean)
•all corrections applied
•comes with error grids
•JPL, CSR and GFZ solutions

### Page 43

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Upcoming mission: GRACE-FO (4/2018)
•FO = “Follow-On”
•NASA/GFZ, launch April 2018
•Like GRACE but in addition carries
technology demonstrator Laser
Range Interference instrument,
•Uncertainty 10 to 20 times  better
than GRACE (theoretically)
•With lessons learnt, somewhat
better resolution than GRACE will be
obtained
•Same data products, same formats
•With some luck we may soon have a
30-year time series of water storage
variability


### Page 44

GES-02 Advanced Data Analysis: Physical Geodesy                                                  2 Intro to GRACE/FO Take-home message
* GRACE provides total (vertically integrated water storage), beginning
  2002 until end of 2017.
* Data is freely available.
* Latency ~ 6 weeks, but some products are available with reduced
   accuracy
* Resolution is ~300 km at monthly scale, around 600 km at daily scale
* Groundwater storage change can be derived from TWS minus
   modeled soil moisture change from hydrological models.
* Higher spatial resolution through assimilation into hydrological models
* Can not measure the total mass  of the ocean or land water
* Can not separate it from the solid Earth mass


---

## lec-02.pdf

_Total pages: 41_

### Page 1

GES-02 Advanced Data Analysis: Physical Geodesy                                                                3 Solid Earth loading: Theory
Makan Karegar
April 28, 2022
karegar@uni-bonn.deSolid Earth Lading: Theory

### Page 2

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory                             Content
* Overview of stress and strain tensors
*  Constitutive equation: relation between strain and stress in an elastic medium
*  Spherical and layered Earth model
*  Spatial convolution of Green’s function and load
*  Boussinesq problem: half-space models
*  Three equations for a spherical Earth model: equation of conservation of
   momentum, Poisson equation and equation of continuity
* Loading Love numbers
* Green’s function
* Spherical harmonic approach for computing displacement

### Page 3

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory                         Stress tensor
We distinguish between body forces  that act at all points in the earth, such as
gravity, from surface forces .
Surface forces are those that act either on actual surfaces within the earth, such
as a fault or an igneous dike, or forces that one part of the earth exerts on an
adjoining part.
The traction (stress) vector  T is defined as the limit of the surface force per unit
area dF acting on a surface element dA, with unit normal n as the size of the area
element tends to zero
                                                         T = lim ∆F/∆A
                                                          ∆A           0


### Page 4

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
The traction depends:

1. Forces acting on the body

2. Orientation of the surface elements
Let’s formulate the traction components acting on three mutually orthogonal surfaces
populate a second-rank stress tensor σ
The mean normal stress  is equal to minus the pressure
σkk/3 = − p    also called volumetric stress .
where
σkk  = σ11 + σ22 +σ33                        Stress tensor

### Page 5

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
Displacement vector:     u = ( u1, u2, u3) or u = (ux, uy, uz)
Cartesian coordinates:  x = ( x1, x2, x3) or x = (x, y, z)
Linearization of the displacement using Taylor’s series expansion as:
we can rewrite the second term of expansion ( .∇u) as:
where we define as symmetric tensor εij as strain tensor and anti symmetric
tensor wij as rotation tensor:
∇.u               Strain and rotation tensor

### Page 6

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
                    Strain tensor
* The strain tensor is symmetric and only 6 components are independent:
   ε21 = ε12 , ε31 = ε13 , ε23 = ε32
* Strains are unitless but have directions,  e.g. axial strain, shear strain...In vector form:

### Page 7

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
Diagonal elements: axial strains Off-diagonal: Shear Strains
Applied stressFor example a dilatation (uniaxial strains ):
i.e. length changes in the x,y or z directions
Strain in x direction
Volume change:  sum
of  strains in x,y,z
directions
=(ε11+ ε22+ ε33)                    Strain tensor
dX/X = εxx dY/Y = εyy

### Page 8

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
Diagonal elements: axial strains Off-diagonal: Shear Strains
Shear strain examples: think of angle change of a side
                    Strain tensor
xx+dx
θdy
shear strain: (dy/x)= tanθ = θ
as displacement in the orthogonal
direction divided by the lengthA simple picture for shear strain

### Page 9

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
General linear elasticity is described by generalized Hooke’s law:
σij    stress tensor
εkl    strain tensor
Cijkl   elastic properties of a material       Constitutive equation: elastic rheology
The constitutive equation describes how material stress and strain (or strain rate) are
related to each other which define a rheology.
On short timescales (seconds to years and several decades) the Earth deforms well
with elasticity but on longer (geologically) large time scales (multi-decades to
thousands and millions of years), the Earth mantle behaves like a very viscous fluid .

### Page 10

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
For an isotropic body  (properties of materials are directional independent) the
elastic parameters reduce to two scalar parameters λ and μ, called Lame’s
parameters . The result is Hooke’s law  (stress and strain or their rates are linearly
dependent).
σij stress tensor       εij  strain tensor
εkk =  (ε11+ ε22+ ε33) or volumetric strain
μ shear modulus (unit: Pa) relating shear stress to shear strain
λ is unitless (no name!) but both called Lame’s parameters
Shear component (i ≠ j) gives:
σij = 2μεij  e.g. i = 1 and j = 2 σ12 = 2με12
Normal components (i = j = 1, 2, or 3):
σii = (3λ + 2μ)εii = 3Kεii       Constitutive equation: elastic rheology

### Page 11

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory Lame’s parameters can define other elastic parameters which are convenient
for particular applications:
Young‘s modulus
(Unit: Pa)
Poisson number
(unitless)
σkk =  (σ11+ σ22+ σ33) or volumetric stress
Bulk or compression modulus ( K) relates mean normal stress σkk/3 to volumetric
strain εkk :
There are five widely used elastic constants ( μ, λ, K, ν, E ) only two are independent.
       Constitutive equation: elastic parameters

### Page 12

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory        Constitutive equation: elastic parameters
                  Some physical insights
Young’s modulus relates stress and strain for the special case of uniaxial
stress:
Strain in direction of applied stress is proportional to 1/ E.
E is always a positive number and it can be viewed as a measure of stiffness!
Applied stress


### Page 13

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory        Constitutive equation: elastic parameters
                  Some physical insights
Poisson’s ratio measures the ratio of strain in the orthogonal direction to that in the
direction the stress is applied.
Applied stress
Strain in direction of normal to stress direction is proportional to  ν
ν is bounded by −1 ≤ ν ≤ 0.5. For ν = 0.5, 1/ K = 0, and the material is incompressible

### Page 14

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
Sandstone:
E = 20 GPa = 20∙109N/m2
Granite:
E = 50 GPa
Strain: q=0.001    Normal stress 5 MPa
   Normal stress 5 MPa
Strain: q=?       Constitutive equation: elastic parameters
                  Some physical insights
The elastic parameters of rocks are computed in labs be applying a known stress and
measuring strain.

### Page 15

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory        Constitutive equation: elastic parameters
1)Measuring travel times of Earthquake-generated seismic waves (P wave and S
wave)
2)Velocities vp and vs can be converted to elastic parameters using empirical
relations such as
Primary (compressional)
 P-WaveSecondary (shear)
S-Wave
This is done in seismology:


### Page 16

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
       Constitutive equation: elastic parameters
                     values within the Earth

### Page 17

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
We want evaluate the deformations at point of interest like at a GPS site with
a coordinate ( θ,λ) at time t due to surface density variations ∆σ(θ’, λ’), that is
at time t over a certain integration domain dσ:
a is average radius of Earth.
Ψ is spherical distance between ( θ,λ) and (θ’, λ’):
G(Ψ) is a Green’s function for vertical or horizontal displacement and is calculated
based on elastic-half space model  or spherical symmetric isotropic layered model .
In general the displacement Green’s function is response of an Elastic earth to
unit point load mass exerted at the pole. load
r
rΨComputation
point (θ,λ)Integration
        point (θ’,λ’)
                                Spatial convolution of Green’s function and load

### Page 18

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory The convolution integral is generally discretized and numerically evaluated
by means of summation of the products  σ(θ’, λ’)G(ψ) over j compartments
within the  entire globe or a specific domain.
(θ’, λ’)(θ’, λ’)(θ, λ)                              Numerical evaluation of convolution integral
Stokes integral is another
convolution integral in geodesy
which is used for geoid
determination (next lectures)

### Page 19

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
  Boussinesq problem
* A cylindrical coordinates system ( r, θ, z),
  centered on the force application point.
* A point load PB acts on the surface
* Stress-free surface
(no shear stress on surface)
* Half-space elastic
(Lame parameters defined the elastic rheology)
* Homogeneous and isotopic Earth
(Elastic parameters are radially and laterally constants)
 μ(r) = μ and  λ(r) = λ
* Non-gravitational Earth (the Earth gravity filed is ignored)
* Flat Earth assumption  Free surface
Boussinesq, J (1885). Applications des Potentials a l’ Equilibre et Mouvement des Solides Elastiques. p.231.

### Page 20

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
* Because the surface point load is axially symmetric, the
  displacement doesn’t depend on the longitude and all
  components of surface displacements are independent of  θ
* There are only two adjustable model parameters:
   the Young’s modulus ( μ) and the Poisson’s ratio ( ν).
* The sensitivity to the Poisson’s ratio is actually negligible so that in practice
only the Young’s modulus is adjusted and ν is generally set to a standard value
of 0.25
* The half-space model may be reasonable if the load extends over small area
   for example lakes, dams, etc.P = Load magnitude for example a column of
water thickness h exerts a pressure P = ρgh
R = horizontal distance between the
       observation point and the point loadDisplacement field in Boussinesq problem
* Analytical  solution

### Page 21

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
ur, uz are Green’s function for vertical and horizontal displacement and should be
convolved with distribution of loads (here surface density) within an area of
interest to calculate the deformation at a point of desire:
Because the extend of loads are small for the half-space model, the convolution
integral is evaluated over a limited domain within which loading data are
available (e.g. within a lakes or a dams and etc.)  Displacement field in Boussinesq problem


### Page 22

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory body Development of half-space models
The half space model can be further developed to take into account the complex
shape of loads and the radial non-homogeneity of elastic parameters:
Some examples:
* A rectangular load  has been solved by Becker & Bevis (2004) and obtained
  a semi-analytical solution.
* A simple line  model by Jaeger et al. (2007) and used first by Jiang et al (2010)
  to model uplift in Greenland due to mass loss.
* A multi-layer half pace model suggested by Pan (1997)
  to take into account the effect of strata


### Page 23

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
For elastic half-space model, load always cause downward and inward
deformationDisplacement field in half-space model


### Page 24

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory    Spherical and layered Earth model
Elastic parameters and density vary only radially and are constant laterally.
This is 1D or spherically symmetric Earth model.
Elastic and Self-gravitating
PREM Earth model (Dziewonski and Anderson, 1981)
Iasp91 Earth model (Kennett and Engdahl, 1991)
ak135 model (Kennett et al., 1995)radius (r)

### Page 25

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
Three equations for a spherical Earth model
1. Equations of motion or equation momentum conservation
2. Poisson equation inside the Earth
3. Continuity equation
After integrating these three differential equations with some boundary condition
for a spherical Earth model, the unitless elastic loading Love numbers  are
computed.
The Green’s functions  are then calculated using these Love numbers, allowing to
calculate the displacement or stress for distributed loads using Green’s integral


### Page 26

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
Equation of conservation of momentum
Body forces F include act at all points in the Earth and include:
1.  Gravity including gravitational & centrifugal forces
2.  Loads including those in physical contact with Earth’s surface such as
     hydrological loads, ice-water redistribution, atmospheric pressure,
     earthquake, etc. or gravitation pull from outside like tidal forces.
σ:  stress tensor
ρ:   density
 Law conservation of momentum  requires that the body forces  F (per unit mass)
acting on the element of the body are balanced by the stresses that act on the
surface of the element .

### Page 27

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
The stress tensor is the sum of the initial hydrostatic pressure  (p0)
and a shear change  (non-hydrostatic stress that later can be related to
strain)
The hydrostatic pressure p0 has negative sign (compressive stress)
* positive means stress act outward normal to the surface.
so equation of conservation of momentum gets:
Let the body deform elastically u = ( uθ, uλ, ur) in t0 , then the pressure
after a small time  increment δt will be:
Note: pressure increase mean negative displacement and
vice versa
then equation of conservation of momentum reads:Equation of conservation of momentum

### Page 28

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory The gradient of initial pressure p0 getsNote the hydrostatic pressure depends only on radius (depth):
[N/m2] or [kg/s2m]
ρ0 Initial density
where er is the unit vector, positive outward from the Earth center.
Equation of conservation of momentum
 p(r) = ρ g r ρ: density
g: gravitational attraction
     both depends on depth


### Page 29

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory Let’s assume that the body force F is only gravity - it can be
expressed as the negative gradient of the potential:
The potential is sum of initial state Φ0 and the small perturbation Φ1
ρF  =   - ρ ∇(Φ0 +Φ1) = - ρ ∇Φ0 -ρ ∇Φ1
The initial density is also changes ρ0 :            ρ = ρ0 +ρ1
ρF = - ρ ∇Φ0 -ρ ∇Φ1 =  - (ρ0 +ρ1) ∇Φ0 - (ρ0 +ρ1) ∇Φ1 =  - ρ0∇Φ0 -ρ1∇Φ0 - ρ0∇Φ1 -ρ1∇Φ1
where from previous slide  ρ0∇Φ0 = ρ0 g   is ∇ρ0   and let  ρ1∇Φ0 = ρ1g    and   ρ1∇Φ1 ≈0
Now replace in momentum equation
and we get:
Equation of conservation of momentum

### Page 30

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory   Contribution due to
  shear stress  changeContribution of
  deformation
Contribution due to changes
in gravitational potential (gravity),
often called self-gravitational effectEffects of density
changes (compressibility )
Equation of conservation of momentum
For a incompressible Earth effects of pressure on the density are zero so the last term
will be zero.
PREM is constructed based on an incompressible model.

### Page 31

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory The perturbed gravitational potential field satisfies the Poisson equation
within the volume of the Earth:
For an incompressible Earth the density change ρ1 will be zero so equation
above reduces to the Laplace equation:
                 Poisson equation

### Page 32

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
Three equations for a spherical Earth model
 Equations of motion:
Poisson equation inside the Earth:
Continuity equation:
For an incompressible Earth model the density perturbation ρ1 will be zero so we
will get two equations:
uSpherical harmonic solution for displacement vector u = urer + uvev and
perturbed potential Φ:
with radial-dependent spherical harmonic coefficients Un (r), Vn(r) and Φn(r).
* Un (r), Vn(r) have unit of meter and Φn(r) has potential unit.


### Page 33

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory Since we want to compute deformation, stress and strain on the Earth’s surface
r = R
The Un (R), Vn(R) and Φn(R) coefficients can be related to potential of point load
mass W and some degree dependent coefficients, the so-called loading Love
numbers as:
                               Un = W hn/g         Vn = W ln/g       Φn = W kn
hn actually converts the potential of point load mass ( W) to coefficients of vertical
displacement
ln actually converts the potential of point load mass ( W) to coefficients of horizontal
displacement
kn actually converts the potential of point load mass ( W) to coefficients of
perturbed potential
g is gravitational acceleration at Earth’s surface                                 Loading Love numbers

### Page 34

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory f
Putting back the Un (R), Vn(R) and Φn(R) coefficients into differential equations and
considering the elastic consecutive equation, with applying the boundary
conditions corresponding to different layers of 1D Earth model, the Love number
are estimated (Farrell, 1972).

hn first Love number ln second Love number kn third Love number
hn is used for calculating
    vertical deformationln is used for calculating
    horizontal deformationkn is used for calculating
        potential (mass)
             changes                       Loading Love numbers

### Page 35

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory                        Use Love numbers to calculate displacement (Green’s
                                                           function)
u Point load
r
rΨUr:vertical disp.Uv: horizontal disp.
For a point load with unit mass:
M mean mass of Earth
R mean radius of Earth

vGPS or
computation point
Similar Green’s functions exist for indirect effect
of gravity, effect of tilt and strain components at
the surface (Farrell, 1972)

### Page 36

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory Use Kummar’s transformation to converge the Green’s functions
v
Convergence problem?
v
v


### Page 37

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory
Loading effect
Third Love number Kl
Total changes in potential is sum of potential load and indirect effect due to solid Earth
deformation. The last term is represented by third Love number Kl.
Φn = W kn
Kl is used to relate the potential and surface density SHCs
Potential (geoid height) SHCs ↔ Mass density SHCs
 What GRACE delivers


### Page 38

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory xSpherical harmonic approach for computing displacement
Note the integration domain is entire
Earth
‘ ‘
By putting Pn in Gu and then in convolution integral and using the spherical
harmonic expression of surface load Δσ:The convolution integral for calculating deformation can be represent by
spherical harmonic expression which is easier to use when dealing with global
loading data such as GRACE and hydrological models:
Potential (geoid) SHCs or often called potential
 Stokes coefficients (what GRACE delivers
to us)

### Page 39

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory after applying orthogonal propertiy of Legendre polynomial, after a few
simplification we get ( see Mitrovia et al. 1994, JGR ): Spherical harmonic approach for computing displacement
λ λ
λ λ
λ λUp component of displacement:
East component of displacement:
North component of displacement:
The expansion begins from degree 1. Degree 0 corresponds to entire mass
change of Earth and is set to zero due to conservation of mass.

### Page 40

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory                             Content
* Overview of stress and strain tensors
*  Constitutive equation: relation between strain and stress in an elastic medium
*  Spherical and layered Earth model
*  Spatial convolution of Green’s function and load
*  Boussinesq problem: half-space models
*  Three equations for a spherical Earth model: equation of conservation of
   momentum, Poisson equation and equation of continuity
* Loading Love numbers
* Green’s function
* Spherical harmonic approach for computing displacement

### Page 41

GES-02 Advanced Data Analysis: Physical Geodesy                                                               3 Solid Earth loading: Theory                             Content
* Overview of stress and strain tensors
*  Constitutive equation: relation between strain and stress in an elastic medium
*  Spherical and layered Earth model
*  Spatial convolution of Green’s function and load
*  Boussinesq problem: half-space models
*  Three equations for a spherical Earth model: equation of conservation of
   momentum, Poisson equation and equation of continuity
* Loading Love numbers
* Green’s function
* Spherical harmonic approach for computing displacement


---

## lec-03.pdf

_Total pages: 31_

### Page 1

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application Advanced Data Analysis
Lecture 5: Solid Earth loading: Applications
Makan Karegar
June 22, 2023
karegar@uni-bonn.de

### Page 2

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application                              Content
Tidal and non-tidal ocean loading
Tidal and non-tidal atmospheric loading

                                                                  Hydrological loading
Cryospheric loading


### Page 3

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Objectives:
1. Reduction (correction)  of geodetic observations for studying a specific process – e.g. correcting GPS time series to find
deformations due to tectonics, volcano, slow-slip ans etc.
2. Understanding  processes (e.g. ocean) from observations, which include the response of the solid Earth (indirect
Effects)
3. Better understanding of the Earth structure , through improving Earth models from observationsLoading: corrections for geodetic data and/or a process to
study Earth’s structure and underlying processes

### Page 4

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application Use of Earth models in geodesy
Elastic Models
- homogeneous solutions of the PDEs >> Love numbers
- Earth tides
- Loading >> mass transports within the Earth system
- Atmospheric loading / gravitation of the atmospheric masses
- Ocean tides / ocean load tides
- non-tidal ocean mass change (e.g. changes in currents, sea level change)
- Snow/ice accumulation and melting / polar regions
- Water storage variability (groundwater, surface water, …)
Anelastic models
- In reality, Earth tides follow the tide-generating body slightly delayed (i.e. with a small lag time or
angle). Another example: Post-seismic deformation
Viscoelastic models (similar to Anelastic models but much longer time scale)
- Postglacial rebound / GIA, Earthquake modeling (post-seismic deformation)
Poroelastic models
- Earthquake geodesy, tectonics, aquifer compaction

### Page 5

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application Forward model:
From mass to displacement
Inverse model:
From displacement to mass
B: vector of GPS observations of the seasonal vertical oscillation, σ: the vector of standard errors,
X: vector of surface water mass at each pixel, A: design matrix consisting of the Green’s functions
L: Laplacian operator, β: a regularization factorForward model vs. Inverse model

### Page 6

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
A simple harmonic model for analyzing loading time series
Linear trendAnnual variation
Semi-annual
variationAnnual variation
Your data
Model parameters are calculated using LSQ.
Alternatively, the RMS scatter of time series around a linear term is use to demonstrate
variability of loading data .

### Page 7

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
                            The ocean tides for harmonic M2
                           OTL provider:  http://holt.oso.chalmers.se/loading/
 Gives you displacement for specific lat/longTidal ocean loading

### Page 8

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Tidal ocean loading

### Page 9

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
RMS of variability of predicted ocean bottom pressure
height displacements for every 2.5 degree
Maximum predicted surface displacement
What caused ocean bottom pressure change?
 1)  Water mass flux: Evaporation and
      Precipitation, Ice melt, river discharge
2)  Ocean circulation
3)  Atmospheric pressure change over ocean
     (inverse barometric (IB) effect)Non-tidal ocean loading or ocean bottom pressure loading

### Page 10

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
van Dam & Ray (2010) S1 and S2 Atmospheric Tide Loading Effects for Geodetic ApplicationsAtmospheric tidal loading
Cosine mode Sine mode
Amplitude
Phase

### Page 11

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Largest changes close to poles
Lowest changes along the equator due to small variability in the solar energyNon-tidal atmospheric loading
2000-2010National Centers for Environmental Prediction (NCEP)
pressure data

### Page 12

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Non-tidal atmospheric loading
* Larger variability at mid-lat to high-lat
* Horizontal displacement is 5 times
smaller than the vertical displacement
*Coastal areas less affected by
atmospheric loading due to the ocean
response that tends to mitigate the
effects of loadingvertical displacement
horizontal displacement

### Page 13

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Where to look for atmospheric and oceanic loading data?
This paper gets you a nice overview!
ESMGFZ: three-hourly 0.5 ◦ ×0.5 ◦ grid of non-
tidal atmospheric and oceanic loading
displacement models
ATMIB: three-hourly 0.5 ◦ ×0.5 ◦ grid of non-tidal
atmospheric loading displacement model
ERAin: six-hourly 0.5 ◦ ×0.5 ◦ non-tidal
atmospheric loading displacement model
ATMMO: six-hourly 0.5 ◦ ×0.5 ◦ tailored non-tidal
atmospheric and oceanic displacement model
ECCO2: daily 0.5 ◦ ×0.5 ◦ non-tidal oceanic
loading displacement model

### Page 14

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Comparing Non-tidal oce and atm models with GPS over Europe

### Page 15

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Hydrological loading
Source: UNAVCOHydrological loading:
* surface water: lakes, reservoirs and rivers
* groundwater
* snow pack and ice
* soil moisture
* vegetationTWS

### Page 16

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Poro-elastic effect:
Charge and recharge of aquifer could
cause poro-elastic deformation: Think of
opposite deformation of elastic loading
Source: UNAVCO

### Page 17

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Half space model: Line model – Greenland uplift


### Page 18

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application River Loading in the Amazon: (Bevis et al. 2005)
Half space model: Rectangular forward model

### Page 19

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Hydrological loading: uplift due to drought in Texas
Karegar et al. (2014)North American Land Data Assimilation System (NLDAS): hourly resolution at 0.125O grid
Data access: https://disc.gsfc.nasa.gov/datasets?keywords=NLDAS

### Page 20

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Half space model: Rectangular inverse model
Uplift at GPS sites 2010-2012Total mass loss due to
droughtRectangular half-space
model
Karegar et al. (2014)

### Page 21

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Half space model: Application to lakes and reservoirs and
probing local structure of crust
The elastic parameters  are
estimated by optimally
fitting the modeled
displacements to observed
displacement (InSAR)
through a half-space model
and a given water load

### Page 22

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Inverting GPS vertical displacement for a spherical layered Earth model
Cosine modeCosine mode Sine modeFrom GPS
From a hydro model
From GRACETWS (m)


### Page 23

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Hydrological loading: GLDAS model
Global Land Data Assimilation System (GLDAS):
Hourly resolution and 0.25O grid
Data access:
https://disc.gsfc.nasa.gov/datasets?keywords=GLDAS

### Page 24

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Hydrological loading: GRACE model
Data access:
Lecture 3

### Page 25

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Hydrological loading: Combining GRACE and a hydro model
Karegar et al. (2018)

### Page 26

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Hydrological loading: Combining GRACE and a hydro model
Karegar et al. (2018)

### Page 27

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Hydrological loading: Assimilating GRACE into a hydro model
Springer et al. (2019)TWS: GRACE TWS: GRACE-assimilated model

### Page 28

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Hydrological loading: Assimilating GRACE into a hydro model
What happens in spectra of loading data:
Springer et al. (2019)

### Page 29

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
Milliner et al. (2018)Elastic response of crust to fast loading: Hurricane Harvey 2017
Water load from rainfall
caused up to 1 cm land subsidence.


### Page 30

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application
 Land subsidence fro GPS is inverted to
estimate the added mass from the rainfallElastic response of crust to fast loading: Hurricane Harvey 2017

### Page 31

GES-02 Advanced Data Analysis: Physical Geodesy                                                       5 Solid Earth loading: Application Earth tides / global (e.g. for satellite orbit computation) elastic Earth
model / Love numbers are sufficient.
Earth tides / local (e.g. station motion, gravimeter) is well-known, local
properties of the Earth model (high SH degrees) become important, could
cause problems for precise measurements.
Loading effects at time scales up to years / decades: elastic Earth models
if sufficient , but loading mass often not well-known!
Loading at longer time scales (postglacial rebound): viscoelastic Earth
rheology not well-known, loads not well-known, this is often a big problem.To what accuracy can we compute loading effects?


---

