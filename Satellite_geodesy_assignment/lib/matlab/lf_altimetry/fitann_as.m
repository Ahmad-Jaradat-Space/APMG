%%Function which fits an annual semi annual and trend and mean to a time
%series 
%input dat is assumed to consist of two or three columns ( unweighted
%and weighted LS resp)
% 

function [fitann]=fitann(dat)

%annual and semi annual frequency (time expressed in years)

%ohm1=2*(365/351.9)*pi;
ohm1=2*pi;
ohm2=4*pi;

ncol=size(dat,2);


ndat=length(dat(:,1));
%subtract approximate center time from time vector
mean=round(sum(dat(:,1))/ndat);
dat(:,1)=dat(:,1)-mean;

%create design matrix

%annual part
A=[cos(ohm1*dat(:,1)) sin(ohm1*dat(:,1))];

%semiannual part (append to right)
A=[A cos(ohm2*dat(:,1)) sin(ohm2*dat(:,1))];

%drift (in m per year) plus mean(append to right)
%A=[A dat(:,1) ones(ndat,1)];

%amount of unknowns
nunk=size(A,2);



%solve
if ncol ==3

P=diag([dat(:,3)].^-2);

N=A'*P*A;

B=A'*P;

else

N=A'*A;
B=A';

end

C=inv(N);

%right hand side of normal system

b=[dat(:,2)];

%left hand solve
SOL=N\(B*b);
%SOL=C*(B*b);

fw=A*SOL;
res=b-fw;

fitann.res=res;
fitann.fw=fw;

%VCE estimation
if ncol==3
  q=(res'*P*res)/(ndat-nunk);
else
  q=(res'*res)/(ndat-nunk);
end

fitann.VCE=q;
SOL=[SOL sqrt(q*diag(C))];



fitann.trend=SOL(nunk-1,1);
fitann.trend_sig=sqrt(SOL(nunk-1,2));
fitann.mean=SOL(nunk,1);
fitann.mean_sig=sqrt(SOL(nunk,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Conversion to plain Cosine harmonic


%calculate cosine amplitude and phase and their errors
% A cos(ohm t)+B sin (ohm *t) =AM cos(ohm*(t-PH))
COSrepr=SOL;
%Annual part
ind1=1;
ind2=2;
Ac=SOL(ind1,1);
As=SOL(ind2,1);
ohm=ohm1;
AM=sqrt(Ac^2+As^2);
PH=atan2(As,Ac)/ohm;

%error propagation
AMsig=Ac^2/(Ac^2+As^2)*diag(C(ind1,ind1))+...
      As^2/(Ac^2+As^2)*diag(C(ind2,ind2))+...
      2*As*Ac/(Ac^2+As^2)*diag(C(ind1,ind2));
AMsig=AMsig*q;

PHsig=(As/((ohm*Ac^2)*(1+(As/Ac)^2)))^2*diag(C(ind1,ind1))+...
      (1/((ohm*Ac)*(1+(As/Ac)^2)))^2*diag(C(ind2,ind2))-...
      2*(As/(Ac^3))*(1/(ohm*(1+(As/Ac)^2)))^2*diag(C(ind1,ind2));

PHsig=PHsig*q;

fitann.an_AM=AM;
fitann.an_AM_sig=sqrt(AMsig);
fitann.an_PH=PH;
fitann.an_PH_sig=sqrt(PHsig);

% fw2(1:ndat,1)=COSrepr(ind(1,1),1)*cos(ohm*(dat(:,1)-COSrepr(ind(2,1),1)));
% fw2(ndat+1:2*ndat,1)=COSrepr(ind(1,2),1)*cos(ohm*(dat(:,1)-COSrepr(ind(2,2),1)));
% fw2(2*ndat+1:3*ndat,1)=COSrepr(ind(1,3),1)*cos(ohm*(dat(:,1)-COSrepr(ind(2,3),1)));


%semi annual part
ind1=3;
ind2=4; 
Ac=SOL(ind1,1);
As=SOL(ind2,1);
ohm=ohm2;
AM=sqrt(Ac^2+As^2);
PH=atan2(As,Ac)/ohm;

%error propagation

AMsig=Ac^2/(Ac^2+As^2)*diag(C(ind1,ind1))+...
      As^2/(Ac^2+As^2)*diag(C(ind2,ind2))+...
      2*As*Ac/(Ac^2+As^2)*diag(C(ind1,ind2));
AMsig=AMsig*q;

PHsig=(As/((ohm*Ac^2)*(1+(As/Ac)^2)))^2*diag(C(ind1,ind1))+...
      (1/((ohm*Ac)*(1+(As/Ac)^2)))^2*diag(C(ind2,ind2))-...
      2*(As/(Ac^3))*(1/(ohm*(1+(As/Ac)^2)))^2*diag(C(ind1,ind2));

PHsig=PHsig*q;


fitann.san_AM=AM;
fitann.san_AM_sig=sqrt(AMsig);
fitann.san_PH=PH;
fitann.san_PH_sig=sqrt(PHsig);




