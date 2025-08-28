function plot_histogram_and_pdf(v, sigma2, component_name)
% PLOT_HISTOGRAM_AND_PDF - Plots residual histogram with overlaid normal PDF
%
% INPUTS:
%   v              - residuals (Nx1)
%   sigma2         - variance (sigma^2), either given or estimated
%   component_name - name of the coordinate component ('East', 'North', 'Up')
%
% OUTPUT:
%   Creates a figure showing the histogram of residuals and a normal distribution overlay.

% Clean residuals
v = v(~isnan(v));
v = v(:); % ensure column vector

if isempty(v)
    warning(['No valid residuals for ' component_name '.']);
    return;
end

% Create figure
figure('Name', [component_name ' - Residual Histogram']);
histogram(v, 'Normalization', 'pdf', 'FaceColor', [0.3 0.3 1]);
hold on;

% Generate normal distribution
x = linspace(min(v), max(v), 200);
pdf_normal = normpdf(x, 0, sqrt(sigma2));
plot(x, pdf_normal, 'r', 'LineWidth', 2);

% Labels and formatting
title([component_name ' Component - Residual Histogram vs Normal PDF']);
xlabel([component_name ' Residuals [mm]']);
ylabel('Probability Density');
legend('Histogram', 'Normal PDF');
grid on;
axis tight;

end
