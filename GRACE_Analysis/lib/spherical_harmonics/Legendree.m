function P = Legendree(lat, L)
c = cosd(lat);
s = sind(lat);
P = zeros(L+1, L+1);
for n = 0:1:L
    for m = 0:1:n
        an = sqrt((2*n + 1) / (2*n));
        bn = sqrt(2*n + 1);
        cn = sqrt((2*n + 1) / ((n - m) * (n + m)));
        dn = sqrt(2*n - 1);
        en = sqrt(((n - m - 1) * (n + m - 1)) / (2*n - 3));
        if (n == 0 && m == 0)
            P(n+1, m+1) = 1;
        elseif (n == 1 && m == 1)
            P(n+1, m+1) = sqrt(3) * s;
        elseif (n == m)
            P(n+1, m+1) = an * s * P(n, m+1);
        elseif (n - m == 1)
            P(n+1, m+1) = bn * c * P(n, m+1);
        else
            P(n+1, m+1) = cn * (dn * c * P(n, m+1) - en * P(n-1, m+1));
        end
    end
end
end