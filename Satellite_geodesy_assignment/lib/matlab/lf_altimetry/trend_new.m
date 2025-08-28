function [T,yfit] = trend_new(ty,x,y,fig);
%Ivan Haigh
%March 2011
%
%Fit a linear or quadratic trend to data
% clear, clc
% s = 6;
% ty = 2;
% load Data_sites;
% x = Y;
% y = DB(:,s);
% i = find(x>=1800 & x<=2009);
% j = find(~isnan(y(i)));
% x = x(i(j));
% y = y(i(j));
% [f] = auto_correlation;
% df = f(s);
% clear DB s  N Y i ns LE L GIA CWR C
% fig = 'Y';
%-------------------------

% Lag-1 Autocorrelation
yd = detrend(y);
ts1 = yd(1:end-1);
ts2 = yd(2:end);
df = abs(corr(ts1,ts2));
%% Calculate trend and standard error
switch ty
    
    %linear trend
    case 1
        X = [x ones(size(y))];
        m = X \ y;
        
        %standard error without auto-correlation
        [n,ncolX] = size(X);
        nu = n-2;                               % Residual degrees of freedom
        yfit = X*m;                             % Predicted responses at each data point.
        r = y-yfit;                             % Residuals.
        rmse = (sqrt(sum(r.^2)))/sqrt(nu);      % Root mean square error.
        mse = rmse^2;                           % Estimator of error variance.
        Xi = inv(X'*X);
        se = [rmse*sqrt(Xi(1,1));rmse*sqrt(Xi(2,2))];
                
        %standard error with auto-correlation
        nu2 = (n-2)*(1-df)/(1+df);                   % Residual degrees of freedom
        rmse2 = (sqrt(sum(r.^2)))/sqrt(nu2);    % Root mean square error.
        mse2 = rmse2^2;                         % Estimator of error variance.
        se2 = [rmse2*sqrt(Xi(1,1));rmse2*sqrt(Xi(2,2))];
        
        %check trend
        [mc,sec,msec] = lscov(X,y);

        %export
        T = [m(1),se(1),se2(1)];        
        
    %Quadractic trend
    case 2
        % trend - Quadtractic
        X = [x.^2 x ones(size(y),1)];
        m = X \ y;
        
        %standard error without auto-correlation
        [n,ncolX] = size(X);
        nu = n-3;                             % Residual degrees of freedom
        yfit = X*m;                           % Predicted responses at each data point.
        r = y-yfit;                           % Residuals.
        rmse = (sqrt(sum(r.^2)))/sqrt(nu);    % Root mean square error.
        mse = rmse^2;                         % Estimator of error variance.
        Xi = inv(X'*X);
        se = [rmse*sqrt(Xi(1,1));rmse*sqrt(Xi(2,2));rmse*sqrt(Xi(3,3))];
        
        %standard error with auto-correlation
        nu2 = (n-2)*(1-df)/(1+df);                  % Residual degrees of freedom
        rmse2 = (sqrt(sum(r.^2)))/sqrt(nu2);    % Root mean square error.
        mse2 = rmse2^2;                         % Estimator of error variance.
        se2 = [rmse2*sqrt(Xi(1,1));rmse2*sqrt(Xi(2,2));;rmse2*sqrt(Xi(3,3))];
        
        %check trend
        [mc,sec,msec] = lscov(X,y);
        
        %export
        T = [m(1)*2,se(1)*2,se2(1)*2];
end
%-------------------------

%% Figure
switch fig
    case 'Y'
        figure('units','normalized','position',[0.1 0.1 0.5 0.8]);
        subplot(211)
        hold on
        plot(x,y);
        plot(x,yfit,'r')
        set(gca,'xlim',[min(x) max(x)])
        title([num2str(T(1),'%6.4f'),' +/- ',num2str(T(2),'%6.4f'),', ',num2str(T(3),'%6.4f')]);
        box on
        
        subplot(212)
        hold on
        plot(x,r,'k');
        box on
        
end
%-------------------------










