function y = pnm(n,teta,dm)

%   P = pnm(n,teta,dm) computes the fully normalized associated 
%   Legendre functions of degree n and order m = 0, 1, ..., n, evaluated
%   for each element of teta.  n must be a scalar integer and  must contain 
%   real values between -90 <= teta <= 90.  
%
%   If teta is a vector, P is an (n+1)-by-L matrix, where L = length(teta).
%   The P(m+1,i) entry corresponds to the fully normalized associated 
%   Legendre function of degree n and order m evaluated at teta(i). 
%
%   In general, the returned array has one more dimension than teta.
%   Each element P(m+1,i,j,k,...) contains the fully normalized associated 
%   Legendre function of degree n and order m evaluated at teta(i,j,k,...).
%
%
%   Note that the first row of P is the Legendre polynomial evaluated at teta 
%   (the m == 0 case).
%
%   Acknowledgment:
%
%   This program is based on a Fortran program by Robert L. Parker,
%   Scripps Institution of Oceanography, Institute for Geophysics and 
%   Planetary Physics, UCSD. February 1993.
%
%
%   Reference:
%     [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%         Functions", Dover Publications, 1965, Ch. 8.
%     [2] J. A. Jacobs, "Geomagnetism", Academic Press, 1987, Ch.4.


if nargin < 3
    error('Not enough input arguments')
elseif nargin > 3
    error('Too many input arguments')    
end

if numel(n) > 1 || ~isreal(n) || n < 0 || n ~= round(n)
    error('N must be a positive scalar integer');
end

x=cosd(teta);
if ~isreal(x) | max(abs(x)) > 1 | ischar(x)
    error('X must be real and in the range (-1,1)')
end

classin = superiorfloat(x);

% The n = 0 case
if n == 0 
    y = ones(size(x),classin);
    return
end

% Convert x to a single row vector
sizex = size(x); 
x = x(:)';

rootn = sqrt(0:2*n);
s = sqrt(1-x.^2);
P = zeros(n+3,length(x),classin);

% Calculate TWOCOT, separating out the x = -1,+1 cases first
twocot = x;

% Evaluate x = +/-1 first to avoid error messages for division by zero
k = find(x==-1);
twocot(k) = Inf(classin);

k = find(x==1);
twocot(k) = -Inf(classin);

k = find(s);
twocot(k) = -2*x(k)./s(k);

% Find values of x,s for which there will be underflow

sn = (-s).^n;
tol = sqrt(realmin(classin));
ind = find(s>0 & abs(sn)<=tol);
if ~isempty(ind)
    % Approx solution of x*ln(x) = y 
    v = 9.2-log(tol)./(n*s(ind));
    w = 1./log(v);
    m1 = 1+n*s(ind).*v.*w.*(1.0058+ w.*(3.819 - w*12.173));
    m1 = min(n, floor(m1));

    % Column-by-column recursion
    for k = 1:length(m1)
        dm1 = m1(k);
        col = ind(k);
        P(dm1:n+1,col) = zeros(size(dm1:n+1))';
        
        % Start recursion with proper sign
        tstart = eps(classin);
        P(dm1,col) = sign(rem(dm1,2)-0.5)*tstart;
        if x(col) < 0
            P(dm1,col) = sign(rem(n+1,2)-0.5)*tstart;
        end
        
         % Recur from m1 to m = 0, accumulating normalizing factor.
        sumsq = tol;
        for m = dm1-2:-1:0
            P(m+1,col) = ((m+1)*twocot(col)*P(m+2,col)- ...
                  rootn(n+m+3)*rootn(n-m)*P(m+3,col))/ ...
                  (rootn(n+m+2)*rootn(n-m+1));
            sumsq = P(m+1,col)^2 + sumsq;
        end
        scale = 1/sqrt(2*sumsq - P(1,col)^2);
        P(1:dm1+1,col) = scale*P(1:dm1+1,col);
    end     % FOR loop
end % small sine IF loop

% Find the values of x,s for which there is no underflow, and for
% which twocot is not infinite (x~=1).

nind = find(x~=1 & abs(sn)>=tol);
if ~isempty(nind)

    % Produce normalization constant for the m = n function
    d = 2:2:2*n;
    c = prod(1-1./d);

    % Use sn = (-s).^n (written above) to write the m = n function
    P(n+1,nind) = sqrt(c)*sn(nind);
    P(n,nind) = P(n+1,nind).*twocot(nind)*n./rootn(end);

    % Recur downwards to m = 0
    for m = n-2:-1:0
        P(m+1,nind) = (P(m+2,nind).*twocot(nind)*(m+1) ...
            -P(m+3,nind)*rootn(n+m+3)*rootn(n-m))/ ...
            (rootn(n+m+2)*rootn(n-m+1));
    end
end % IF loop

y = P(1:n+1,:);

% Polar argument   (x = +-1)
s0 = find(s == 0);
y(1,s0) = x(s0).^n;

% Calculate the fully-normalized functions.
% For m = 1,...,n, multiply by sqrt(2*(2*n+1))
% And for m = 0, multiply by sqrt(2*n+1)
y(2:n+1,:) = sqrt(2*(2*n+1))*y(2:n+1,:);
y(1,:) = sqrt(2*n+1)*y(1,:);
% For m = 1,3,5,...,n, multiply by -1
y(2:2:n+1,:) = -y(2:2:n+1,:);

% Restore original dimensions.
if length(sizex) > 2 || min(sizex) > 1
    y = reshape(y,[n+1 sizex]);
end

if dm==0
    
elseif dm==1
y(1,:)=[];
y(end+1,:)=0;

elseif dm==2
y(1:2,:)=[];
y(end+1:end+3,:)=0;

elseif dm==-1
y=[y(2,:).*(-1);y];

elseif dm==-2 && n==1
y=[zeros(1,max(size(teta)));y(2,:).*(-1);y(1,:)];

elseif dm==-2 && n >1
y=[y(3,:);y(2,:).*(-1);y(1:end-1,:)];

else 
    error('Incorrect dm')
end

