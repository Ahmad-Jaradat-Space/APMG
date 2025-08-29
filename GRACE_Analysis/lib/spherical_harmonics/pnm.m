function P = pnm(lat, L)
%PNM Fully normalized associated Legendre functions P̄nm at latitude (degrees)
% Input:
%   lat - geodetic latitude in degrees
%   L   - maximum degree
% Output:
%   P   - (L+1) x (L+1) matrix, P(n+1,m+1) = P̄nm(lat)

c = cosd(lat);
s = sind(lat);
P = zeros(L+1, L+1);
for n = 0:L
    for m = 0:n
        an = sqrt((2*n + 1) / (2*n));
        bn = sqrt(2*n + 1);
        if n > 0 && (n - m) > 0
            cn = sqrt((2*n + 1) / ((n - m) * (n + m)));
            dn = sqrt(2*n - 1);
        else
            cn = 0;
            dn = 0;
        end
        if n > 1 && (n - m) > 1
            en = sqrt(((n - m - 1) * (n + m - 1)) / (2*n - 3));
        else
            en = 0;
        end
        if (n == 0 && m == 0)
            P(n+1, m+1) = 1;
        elseif (n == 1 && m == 1)
            P(n+1, m+1) = sqrt(3) * s;
        elseif (n == m)
            % Use previous degree same order (n-1,m-1)
            P(n+1, m+1) = an * s * P(n, m);
        elseif (n - m == 1)
            % Use previous degree same order (n-1,m)
            P(n+1, m+1) = bn * c * P(n, m+1);
        else
            % Two-step recurrence
            P(n+1, m+1) = cn * (dn * c * P(n, m+1) - en * P(n-1, m+1));
        end
    end
end
end