function [h_n, l_n, k_n] = loadLoveNumbers(nmax)

S = load('data/love_numbers/PREM_Love_Numbers_n60.mat');

h_n = S.h_n(1:nmax+1);
l_n = S.l_n(1:nmax+1);
k_n = S.k_n(1:nmax+1);

% Detect whether k_n is tidal (typically positive around degree 2) and convert to load
% For load Love numbers, k'_n (n>=2) is typically negative (~ -0.3 for n=2).
% Assume tidal k; convert to load: k' = k - h
k_n = k_n - h_n;


end