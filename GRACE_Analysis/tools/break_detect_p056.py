#!/usr/bin/env python3
import os
import sys
import math

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fn = os.path.join(ROOT, 'data', 'gps', 'P056.tenv3')
if not os.path.exists(fn):
    print('File not found:', fn)
    sys.exit(1)

t_year = []
y_m = []
with open(fn, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('site'):
            continue
        parts = line.split()
        if len(parts) < 13:
            continue
        try:
            t = float(parts[2])
            up = float(parts[12])  # meters
        except ValueError:
            continue
        if math.isfinite(t) and math.isfinite(up):
            t_year.append(t)
            y_m.append(up)

import numpy as np

if len(y_m) < 100:
    print('Insufficient data points:', len(y_m))
    sys.exit(0)

t = np.array(t_year)
y = np.array(y_m) * 1000.0  # convert to mm
N = len(y)

# Single linear fit
A0 = np.vstack([np.ones(N), t]).T
coeff0, *_ = np.linalg.lstsq(A0, y, rcond=None)
res0 = y - A0 @ coeff0
SSE0 = float(np.sum(res0**2))

# Two-segment search
min_seg = max(30, int(0.1*N))
idx_min = min_seg
idx_max = N - min_seg
SSE_best = float('inf')
best_k = None
c1_best = None
c2_best = None

for k in range(idx_min, idx_max):
    t1 = t[:k]; y1 = y[:k]
    t2 = t[k:]; y2 = y[k:]
    A1 = np.vstack([np.ones(len(t1)), t1]).T
    A2 = np.vstack([np.ones(len(t2)), t2]).T
    c1, *_ = np.linalg.lstsq(A1, y1, rcond=None)
    c2, *_ = np.linalg.lstsq(A2, y2, rcond=None)
    SSE = float(np.sum((y1 - A1 @ c1)**2) + np.sum((y2 - A2 @ c2)**2))
    if SSE < SSE_best:
        SSE_best = SSE
        best_k = k
        c1_best = c1
        c2_best = c2

# Chow F-test
p_restr = 2
p_unrestr = 4
num = (SSE0 - SSE_best) / (p_unrestr - p_restr)
den = SSE_best / (N - p_unrestr)
F = num / den if den > 0 else float('inf')

# p-value via SciPy if available; else approximate using survival function for F
pval = None
try:
    from mpmath import quad
    from mpmath import beta as mpbeta
    d1 = p_unrestr - p_restr
    d2 = N - p_unrestr
    # F-distribution survival function using incomplete beta
    x = (d1*F) / (d1*F + d2)
    # regularized incomplete beta I_{1-x}(d2/2, d1/2) = P(F > F0)
    # Using mpmath betainc
    from mpmath import betainc
    pval = float(betainc(d2/2, d1/2, 0, 1-x, regularized=True))
except Exception:
    pval = None

print('Station P056 structural break analysis')
print(f'  N = {N} points')
print(f'  Best break index: {best_k} / {N}')
print(f'  Best break time (decimal year): {t[best_k]:.4f}')
print(f'  Slope before break:  {c1_best[1]:.3f} mm/yr')
print(f'  Slope after break:   {c2_best[1]:.3f} mm/yr')
print(f'  SSE single: {SSE0:.3e}, SSE two-segment: {SSE_best:.3e}')
print(f'  F = {F:.2f}')
if pval is not None:
    print(f'  p-value (approx): {pval:.3g}')
else:
    print('  p-value: (scipy/mpmath not available) use F comparison to critical values')

if F > 10:  # rough significance threshold
    print('Conclusion: Significant break detected. Two linear trends with an offset are appropriate.')
else:
    print('Conclusion: Break not strongly supported at this threshold.') 