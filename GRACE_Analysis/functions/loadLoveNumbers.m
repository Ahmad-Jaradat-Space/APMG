function [h_n, l_n, k_n] = loadLoveNumbers(nmax, model)
h_n = zeros(nmax + 1, 1);
l_n = zeros(nmax + 1, 1);
k_n = zeros(nmax + 1, 1);
h_n(1) = 0;
l_n(1) = 0;
k_n(1) = 0;
if nmax >= 1
    h_n(2) = 0;
    l_n(2) = 0;
    k_n(2) = 0.021;
end
for n = 2:nmax
    if n == 2
        h_n(n+1) = 0.6149;
        l_n(n+1) = 0.0839;
        k_n(n+1) = 0.3020;
    else
        h_n(n+1) = 0.6149 * (1 + 5/(2*n+1) - 7/(2*n+1)^2);
        l_n(n+1) = 0.0839 * (1 + 3/(2*n+1) - 5/(2*n+1)^2);
        k_n(n+1) = 0.3020 * (1 + 2/(2*n+1) - 3/(2*n+1)^2);
    end
end
end