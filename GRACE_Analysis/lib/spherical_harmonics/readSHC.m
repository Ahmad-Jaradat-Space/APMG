function [cnm, snm] = readSHC(file)
n = file(:, 1);
m = file(:, 2);
nmax = max(n);
cnm = zeros(nmax + 1, nmax + 1);
snm = zeros(nmax + 1, nmax + 1);
for i = 1:size(file, 1)
    degree = n(i) + 1;
    order = m(i) + 1;
    cnm(degree, order) = file(i, 3);
    snm(degree, order) = file(i, 4);
end
end