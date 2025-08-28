%% PROJECT- Modeling the elastic response of Earth’s crust to hydrological loading

% clear command window, variables, close all windows, stop in debugging mode in case of errors
clc; clearvars; close all; dbstop if error;

% set some default values for plot formatting
set(groot,'DefaultTextFontSize', 20,...
    'DefaultAxesFontSize', 20,...
    'DefaultAxesTitleFontWeight', 'bold',...
    'DefaultAxesTitleFontSizeMultiplier', 1,...
    'DefaultAxesXMinorTick', 'on',...
    'DefaultAxesYMinorTick', 'on',...
    'DefaultAxesZMinorTick', 'on',...
    'DefaultTextFontName', 'Arial',...
    'DefaultLineLineWidth',3,...
    'DefaultLineMarkerSize',7,...
    'DefaultStemLineWidth',3,...
    'DefaultStemMarkerSize',7,...
    'DefaultTextInterpreter', 'tex', ...
    'DefaultSurfaceEdgeColor', 'flat');

% add all subfolders of the current directory to the PATH
addpath(genpath(pwd));

%% Load data
data = load_tenv3('P140.tenv3');

% extract vectors from the 'data' struct
t = data.t;
% time vector for plotting
dt = datetime(t, 'ConvertFrom', 'modifiedjuliandate');

% number of observations
n = length(t);

% plot of raw data
fig1 = figure('Name', 'P - raw data');
subplot(3,1,1);
plot(dt, data.east, 'b');
title('East Component');
grid on;
xlabel('Time');
ylabel('Displacement [m]');

subplot(3,1,2);
plot(dt, data.north, 'g');
title('North Component');
grid on;
xlabel('Time');
ylabel('Displacement [m]');

subplot(3,1,3);
plot(dt, data.up, 'r');
title('Up Component');
grid on;
xlabel('Time');
ylabel('Displacement [m]');


%% 1.1 estimation of a linear trend function

% max degree of polynomial
p = 2;

% given variance
sigma0_up = mean(data.sig_um);
sigma0_east = mean(data.sig_em);
sigma0_north = mean(data.sig_nm);

% covariance matrices of the raw observations 
Cup = sigma0_up^2 * speye(length(t));
Ceast = sigma0_east^2 * speye(length(t));
Cnorth = sigma0_north^2 * speye(length(t));

% fit of polynomial degree 2
up = fitPolynomial(t, data.up, p, sigma0_up);
east = fitPolynomial(t, data.east, p, sigma0_east);
north = fitPolynomial(t, data.north, p, sigma0_north);

%% 1.2 Plot of the estimated trend functions
figure(fig1);

subplot(3,1,1);
hold on;
plot(dt, east.lest, 'k--');
legend('Raw', 'Fitted');

subplot(3,1,2);
hold on;
plot(dt, north.lest, 'k--');
legend('Raw', 'Fitted');

subplot(3,1,3);
hold on;
plot(dt, up.lest, 'k--');
legend('Raw', 'Fitted');