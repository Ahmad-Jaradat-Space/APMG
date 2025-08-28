function u_vertical = graceToVerticalDeformation(cnm, snm, theta, lambda, h_n, k_n)
constants = physicalConstants();
nmax = size(cnm, 1) - 1;
[nlat, nlon] = size(theta);
u_vertical = zeros(nlat, nlon);
lambda_vec = lambda(1, :);
cosm = cos([0:nmax]' * lambda_vec);
sinm = sin([0:nmax]' * lambda_vec);
for i = 1:nlat
    current_theta_deg = theta(i, 1) * 180/pi;
    Pnm_matrix = Legendree(current_theta_deg, nmax);
    deformation_lat = zeros(1, nlon);
    for n = 1:nmax
        love_factor = h_n(n+1) / (1 + k_n(n+1));
        for m = 0:n
            c_nm = cnm(n+1, m+1);
            s_nm = snm(n+1, m+1);
            Pnm_val = Pnm_matrix(n+1, m+1);
            if m == 0
                trig_terms = c_nm * ones(1, nlon);
            else
                cos_terms = c_nm * cosm(m+1, :);
                sin_terms = s_nm * sinm(m+1, :);
                trig_terms = cos_terms + sin_terms;
            end
            contribution = love_factor * Pnm_val * trig_terms;
            deformation_lat = deformation_lat + contribution;
        end
    end
    u_vertical(i, :) = deformation_lat;
end
u_vertical = constants.R * (constants.rho_water / constants.rho_earth) * u_vertical;
end