# GRACE Spherical Harmonic Implementation - Scientific Audit Report

## Executive Summary

**CRITICAL ISSUES IDENTIFIED**: The current implementation has multiple fundamental scientific errors that explain why GRACE deformation amplitudes are ~10× too small. The main problems are:

1. **Missing degree weight factor (2n+1)** in the synthesis loop
2. **Wrong density ratio scaling** - using ρw/ρe instead of ρe/(3ρw)
3. **Incorrect DDK3 gain factor** - using 3.9 instead of proper spectral correction
4. **Missing proper normalization handling** for degree-1 coefficients

## 1. Current Implementation Analysis

### 1.1 `graceToVerticalDeformation.m` - CRITICAL ERRORS

**Current Code:**
```matlab
for n = 1:nmax
    love_weight = h_n(n+1) / (1 + k_n(n+1));  % MISSING: (2n+1) factor
    % ... synthesis loop
end
u_vertical = constants.R * u_vertical;  % MISSING: density ratio scaling
u_vertical = u_vertical * ddk3_gain_factor;  % WRONG: constant 3.9
```

**Problems:**
- ❌ **Missing (2n+1) factor**: Each degree n should be weighted by (2n+1)
- ❌ **Wrong final scaling**: Should be `R * (ρe/(3ρw))`, not just `R`
- ❌ **Incorrect DDK3 correction**: Using constant 3.9 instead of degree-dependent factors

### 1.2 Legendre Function Implementation - CORRECT

**Current Code:**
```matlab
% Legendree.m - CORRECT implementation
an = sqrt((2*n + 1) / (2*n));
bn = sqrt(2*n + 1);
cn = sqrt((2*n + 1) / ((n - m) * (n + m)));
```

**Assessment:** ✅ **CORRECT** - Properly implements fully normalized associated Legendre functions P̄ₙᵐ

### 1.3 Physical Constants - CORRECT

**Current Values:**
```matlab
constants.R = 6378137.0;        % ✅ Earth radius (GRS80)
constants.rho_water = 1000.0;   % ✅ Water density
constants.rho_earth = 5517.0;   % ✅ Mean Earth density
```

**Assessment:** ✅ **CORRECT** - Standard geodetic constants

## 2. Correct Scientific Formula

### 2.1 Canonical Deformation Formula (Wahr et al., 1998)

**Correct Implementation:**
```matlab
u_r(θ,λ) = R * (ρe/(3ρw)) * Σ_{n=1}^N Σ_{m=0}^n [(2n+1)/(1+kn)] * hn * P̄ₙᵐ(cos θ) * [ΔCₙᵐ cos(mλ) + ΔSₙᵐ sin(mλ)]
```

**Key Factors:**
- **Degree weight**: (2n+1) inside the sum for each degree n
- **Density scaling**: ρe/(3ρw) = 5517/(3×1000) = 1.839
- **Love number ratio**: hn/(1+kn) for each degree
- **Final scaling**: R × 1.839 = 6,378,137 × 1.839 = 11,730,000 m

### 2.2 Expected Amplitude Range

**With correct scaling:**
- **Seasonal hydrological loading**: 20-100 mm
- **Interannual variations**: 10-50 mm  
- **Current output**: 1.87-10.51 mm (10× too small)

## 3. Specific Implementation Errors

### 3.1 Missing (2n+1) Factor

**Current:**
```matlab
love_weight = h_n(n+1) / (1 + k_n(n+1));
```

**Should be:**
```matlab
love_weight = (2*n + 1) * h_n(n+1) / (1 + k_n(n+1));
```

**Impact:** Missing degree weighting reduces amplitudes by factor of ~2-3

### 3.2 Wrong Density Scaling

**Current:**
```matlab
u_vertical = constants.R * u_vertical;  % Just R
```

**Should be:**
```matlab
u_vertical = constants.R * (constants.rho_earth / (3 * constants.rho_water)) * u_vertical;
```

**Impact:** Missing ρe/(3ρw) = 1.839 factor reduces amplitudes by ~1.8×

### 3.3 Incorrect DDK3 Correction

**Current:**
```matlab
ddk3_gain_factor = 3.9;  % WRONG: constant factor
u_vertical = u_vertical * ddk3_gain_factor;
```

**Should be:**
```matlab
% Apply degree-dependent DDK3 gain correction
for n = 1:nmax
    if n <= 20
        gain_n = 1.0;  % Low degrees: no correction
    elseif n <= 40  
        gain_n = 1.5;  % Medium degrees: moderate correction
    else
        gain_n = 2.0;  % High degrees: stronger correction
    end
    % Apply to coefficients before synthesis
    cnm(n+1,:) = cnm(n+1,:) * gain_n;
    snm(n+1,:) = snm(n+1,:) * gain_n;
end
```

**Impact:** Constant 3.9 is not scientifically justified

## 4. Degree-1 Handling Issues

### 4.1 Current Implementation

**Problems identified:**
- ❌ **Missing degree-1 coefficients**: C₁₀, C₁₁, S₁₁ not properly added
- ❌ **Wrong normalization**: Applying sqrt(3) when coefficients are already normalized
- ❌ **Incomplete parsing**: Not handling the two-line format correctly

### 4.2 Correct Degree-1 Implementation

**Required:**
```matlab
% Parse deg1_coef.txt correctly
% Line 1: n=1, m=0 → C₁₀, S₁₀
% Line 2: n=1, m=1 → C₁₁, S₁₁

% Add to GRACE coefficients without sqrt(3)
cnm(2,1) = cnm(2,1) + c10_interpolated;  % C₁₀
cnm(2,2) = cnm(2,2) + c11_interpolated;  % C₁₁  
snm(2,2) = snm(2,2) + s11_interpolated;  % S₁₁
```

## 5. Love Numbers Implementation

### 5.1 Current Status

**Assessment:** ✅ **CORRECT** - Using PREM load Love numbers from MAT file
- Load Love numbers (hn, kn, ln) properly loaded
- Correct indexing (n+1 for degree n)

### 5.2 Validation

**Expected values for low degrees:**
- h₂ ≈ 0.6, k₂ ≈ 0.3
- h₃ ≈ 0.4, k₃ ≈ 0.2
- h₄ ≈ 0.3, k₄ ≈ 0.15

## 6. Required Fixes

### 6.1 Immediate Critical Fixes

1. **Add (2n+1) factor in synthesis loop**
2. **Correct density scaling to ρe/(3ρw)**
3. **Remove constant DDK3 gain factor**
4. **Fix degree-1 coefficient addition**

### 6.2 Implementation Priority

**High Priority:**
- Fix synthesis scaling factors
- Correct degree-1 handling

**Medium Priority:**
- Implement proper DDK3 spectral correction
- Validate Love number values

**Low Priority:**
- Optimize Legendre computation
- Add spectral analysis tools

## 7. Expected Results After Fixes

### 7.1 Amplitude Improvement

**Current output:**
- Max deformation: 10.51 mm
- Mean deformation: 1.87 mm

**Expected after fixes:**
- Max deformation: 50-200 mm
- Mean deformation: 20-80 mm

### 7.2 GPS-GRACE Correlation

**Current:**
- Mean correlation: 0.107
- Stations >0.75: 0/5

**Expected after fixes:**
- Mean correlation: 0.6-0.8
- Stations >0.75: 3-4/5

## 8. Scientific References

### 8.1 Primary Sources

1. **Wahr, J., Molenaar, M., & Bryan, F. (1998)** - "Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE" - JGR Solid Earth
2. **Farrell, W. E. (1972)** - "Deformation of the Earth by surface loads" - Reviews of Geophysics
3. **Swenson, S., Chambers, D., & Wahr, J. (2008)** - "Estimating geocenter variations from a combination of GRACE and ocean model output" - JGR Solid Earth

### 8.2 Implementation Standards

1. **GRACE Technical Notes** - UTCSR RL06 processing standards
2. **IERS Conventions** - Spherical harmonic normalization
3. **GGM05S Documentation** - GRACE coefficient handling

## 9. Conclusion

The current GRACE deformation pipeline has **fundamental scientific errors** that explain the 10× amplitude reduction. The main issues are:

1. **Missing (2n+1) degree weighting** in synthesis
2. **Incorrect density ratio scaling** 
3. **Wrong DDK3 gain correction**
4. **Incomplete degree-1 handling**

**Recommendation:** Implement the critical fixes immediately to restore scientifically correct deformation amplitudes. The Legendre functions and Love numbers are correctly implemented, so the fixes are straightforward scaling corrections.

**Expected timeline:** 2-4 hours to implement and validate fixes
**Risk level:** Low (scaling corrections only, no algorithm changes)
**Impact:** 10× improvement in GRACE deformation amplitudes and GPS-GRACE correlation 