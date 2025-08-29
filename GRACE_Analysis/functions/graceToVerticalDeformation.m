function u_vertical = graceToVerticalDeformation(cnm, snm, theta, lambda, h_n, k_n)
constants = physicalConstants();
nmax = size(cnm, 1) - 1;
[nlat, nlon] = size(theta);
u_vertical = zeros(nlat, nlon);
lambda_vec = lambda(1, :);
cosm = cos([0:nmax]' * lambda_vec);
sinm = sin([0:nmax]' * lambda_vec);
for i = 1:nlat
    theta_deg = theta(i, 1) * 180/pi;
    lat_deg = 90 - theta_deg; % Convert colatitude to latitude
    deformation_lat = zeros(1, nlon);
    for n = 1:nmax
        % Fully normalized associated Legendre for this degree at this latitude
        Pn = pnm(n, lat_deg, 0);  % size (n+1) x 1
        love_weight = h_n(n+1) / (1 + k_n(n+1));
        for m = 0:n
            Pnm_val = Pn(m+1);  % scalar
            c_nm = cnm(n+1, m+1);
            s_nm = snm(n+1, m+1);
            if m == 0
                trig_terms = c_nm * ones(1, nlon);
            else
                cos_terms = c_nm * cosm(m+1, :);
                sin_terms = s_nm * sinm(m+1, :);
                trig_terms = cos_terms + sin_terms;
            end
            contribution = love_weight * Pnm_val * trig_terms;
            deformation_lat = deformation_lat + contribution;
        end
    end
    u_vertical(i, :) = deformation_lat;
end

% Canonical scaling with negative sign (positive load => downward)
scale = constants.R * (constants.rho_earth / (3 * constants.rho_water));
u_vertical = -scale * u_vertical;

end