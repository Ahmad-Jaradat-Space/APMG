## GRACE Deformation Pipeline – Scientific Fix Plan

### Objective
- Restore scientifically correct vertical deformation amplitudes from GRACE by fixing scaling, Love numbers, degree‑1 handling, and baseline, then validate against GPS.

### Diagnosed Issues (What’s wrong)
- **Scaling in synthesis (dominant error):**
  - Missing degree weight (2n+1)/3 and wrong density ratio. Current code multiplies the summed field by `R * (ρw/ρe)` and uses only `h_n/(1+k_n)` inside the sum.
  - Correct physics (fully normalized harmonics):
    - u(θ, λ) = a · [ρe/(3ρw)] · Σ_{n=1..N} Σ_{m=0..n} h′ₙ · [(2n+1)/(1+k′ₙ)] · P̄ₙᵐ(cos θ) · [C̄ₙᵐ cos mλ + S̄ₙᵐ sin mλ].
- **Degree‑1 handling:**
  - `deg1_coef.txt` values are already GRACE‑normalized per line (n=1,m=0 and n=1,m=1 separately). The current parser conflates lines and applies an extra `sqrt(3)` (incorrect), and sets C11=0.
- **Love numbers:**
  - `loadLoveNumbers.m` uses ad‑hoc sequences; must be replaced with real PREM load Love numbers (h′ₙ, k′ₙ, l′ₙ) from `GRACE_Analysis/data/love_numbers/PREM_Love_Numbers_n60.mat`.
- **Baseline and ad‑hoc amplification:**
  - 3‑month baseline can suppress signal; `signal_amplification = 20` is not a scientific deconvolution of DDK3 attenuation.
- **C20 replacement & normalization:**
  - C20 usage is consistent (normalized); keep but re‑validate against time tags.
- **Legendre normalization:**
  - `Legendree.m` computes fully normalized P̄ₙᵐ; keep as‑is.

### Changes to Implement

#### 1) Fix scaling in `GRACE_Analysis/functions/graceToVerticalDeformation.m`
- Inside the n‑loop, include the degree weight (2n+1). Keep h′ₙ/(1+k′ₙ) but multiply by (2n+1):
  - Replace per‑degree factor with: `love_weight = (2*n + 1) * h_n(n+1) / (1 + k_n(n+1));`
  - Use `love_weight` in the contribution accumulation.
- Replace the final scale factor:
  - Remove: `R * (ρw/ρe)`
  - Add:    `R * (ρe/(3*ρw))`
- Ensure the sum starts at n=1 (skip n=0) as currently implemented.

#### 2) Use real PREM load Love numbers in `loadLoveNumbers.m`
- Replace the ad‑hoc generator with a loader for PREM arrays:
  - Load from `data/love_numbers/PREM_Love_Numbers_n60.mat`.
  - Map variables to `h_n`, `k_n`, `l_n` (adjust variable names to the .mat contents; commonly `h`, `k`, `l`).
  - Ensure vectors are length ≥ nmax+1 and indexed so that index (n+1) corresponds to degree n.
  - Return load Love numbers (not tidal), i.e., h′ₙ, k′ₙ, l′ₙ.

Example sketch (adjust variable names to the MAT file):
```matlab
function [h_n, l_n, k_n] = loadLoveNumbers(nmax, model)
  data = load(fullfile('data','love_numbers','PREM_Love_Numbers_n60.mat'));
  % Map variable names to h_n, k_n, l_n as present in MAT
  h_all = data.h; k_all = data.k; l_all = data.l; % verify
  max_available = min([length(h_all), length(k_all), length(l_all)]);
  N = min(nmax+1, max_available);
  h_n = h_all(1:N); k_n = k_all(1:N); l_n = l_all(1:N);
  % Pad if needed
  if length(h_n) < nmax+1
    h_n(end+1:nmax+1) = h_n(end);
    k_n(end+1:nmax+1) = k_n(end);
    l_n(end+1:nmax+1) = l_n(end);
  end
end
```

#### 3) Correct degree‑1 ingestion in `GRACE_Analysis/functions/processGRACEfiles.m`
- Parse `deg1_coef.txt` by filtering on n and m per line (the file contains two lines per month for degree‑1):
  - For lines with n=1, m=0 → C̄10, S̄10 (S̄10 is typically ~0 in this product).
  - For lines with n=1, m=1 → C̄11, S̄11.
- Build time series keyed by month (or decimal year) with all three: C̄10, C̄11, S̄11.
- Interpolate each coefficient to the GRACE mid‑times.
- Assign to SHCs without any `sqrt(3)` (already normalized):
  - `cnm(2,1) = C̄10;`
  - `cnm(2,2) = C̄11;`
  - `snm(2,2) = S̄11;`
- Remove the lines that set C11=0 and the `* sqrt(3)` multipliers.

Parsing sketch:
```matlab
% After reading numeric lines from deg1_coef.txt into rows [YYYYMM, n, m, C, S, ...]
mask10 = rows(:,2)==1 & rows(:,3)==0; % n=1,m=0
mask11 = rows(:,2)==1 & rows(:,3)==1; % n=1,m=1
[t10, c10, s10] = deal(rows(mask10,1), rows(mask10,4), rows(mask10,5));
[t11, c11, s11] = deal(rows(mask11,1), rows(mask11,4), rows(mask11,5));
% Convert YYYYMM to decimal year, then interp c10,c11,s11 to GRACE times.
```

#### 4) Baseline handling in `GRACE_Analysis/main_grace_analysis.m`
- Replace the 3‑month mean with a long‑term static field (e.g., 2004–2009):
  - Use existing `computeStaticReferenceField` with `reference_years = 2004:2009`.
  - For each month, compute `ΔC = C − C_static`, `ΔS = S − S_static`.
- Remove the ad‑hoc amplification:
  - Delete `signal_amplification = 20` and the multiplication of the field by this factor.
- Optional alternative: remove the full temporal mean across the available time span to isolate variability, if a static field isn’t precomputed.

#### 5) Optional: DDK3 attenuation correction (spectral)
- If needed, apply published degree‑dependent gain factors Gₙ for DDK3 in spectral space prior to synthesis:
  - For each degree n: scale ΔC̄ₙᵐ, ΔS̄ₙᵐ by 1/Gₙ (where stable); avoid amplifying noisy high‑n.
  - This replaces any constant amplification and preserves phase.

### Validation & QA
- **Unit sanity:** Deformation returned by synthesis in meters; plots convert to mm.
- **Test point check:** Use `debug_raw_grace.m` at a hydrologically active site; expect seasonal amplitudes of a few mm to ~10 mm after corrections.
- **GPS vs GRACE overlay:** Compare at stations; amplitude ratio ≈ 0.5–1.0 depending on site; correlation improved; reduced bias.
- **Low‑degree signal:** Verify degree‑1 time series are non‑zero and smooth; C20 replacement time series consistent.
- **Spectral check:** Optional—plot average power vs degree pre/post DDK gain correction.

### File Touch List
- `GRACE_Analysis/functions/graceToVerticalDeformation.m` (scaling and degree weight)
- `GRACE_Analysis/functions/loadLoveNumbers.m` (load PREM h′ₙ, k′ₙ, l′ₙ from MAT)
- `GRACE_Analysis/functions/processGRACEfiles.m` (degree‑1 parsing; remove sqrt(3))
- `GRACE_Analysis/main_grace_analysis.m` (baseline via static field; remove amplification)
- `GRACE_Analysis/functions/computeStaticReferenceField.m` (validate normalization/years)

### Rollout Plan
1. Implement scaling fix in `graceToVerticalDeformation.m`.
2. Swap in PREM load Love numbers loader; verify shapes to nmax.
3. Rework degree‑1 parsing and assignment; remove sqrt(3).
4. Change baseline strategy; remove amplification.
5. Quick test at a known site; confirm mm‑level seasonal amplitudes.
6. Full rerun; compare GPS vs GRACE plots; compute amplitude ratios and correlations.
7. If needed, add DDK3 gain correction and re‑validate.

### Rollback
- Revert individual files if unexpected behavior; no schema changes.

### References (for implementation fidelity)
- Wahr, M., Molenaar, M., & Bryan, F. (1998) – loading and GRACE scaling.
- Swenson, S., Chambers, D., & Wahr, J. (2008) – degree‑1 estimation.
- Kusche, J. et al. (2009) – DDK filters and attenuation behavior. 