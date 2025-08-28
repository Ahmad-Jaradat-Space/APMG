clc;
clear all ;
format long;
%% constants
GM = 3.986e14;
R = 6378137;
%% reading the files
addpath("data/Compulsory_gravity/")
GRS80 = load('GRS80.txt') ;
ITSG = load('ITSG_Grace2018.txt');
asin = load('asin.txt');
acos = load('acos.txt');
coast= load('coast.txt');

[ n, cnmg , snmg ]  = readSHC (GRS80);
[ n,cnm , snm ]   = readSHC (ITSG);
[n,cnmsin,snmsin] = readSHC (asin);
[n,cnmcos,snmcos] = readSHC (acos);
n_itsg =   size(cnm,1)-1  ; % max degree for itsg
n_grs =   size(cnmg,1)-1  ; % max degree for grs80
n_sin =   size(cnmsin,1)-1; % max degree for sin file
n_cos =   size(cnmcos,1)-1; % max degree for cos file

%the reslution of the grid = 1 degree
res = 361;
lambda = linspace(-180,180,res);
theta = linspace(-90,90,ceil(res/2));
[L,T]= meshgrid(lambda,theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Task 7 
%Estimate the new SHCs 
% 
%given epochs
t=[0,90,180,273]; % [01/01, 01/04,01/07,01/10]/2010
t0 = 0;

[ ecnm1 ,esnm1] = estimate(cnm , snm , cnmcos , cnmsin , snmcos , snmsin , t0 ,t(1),n_cos );
[ ecnm2 ,esnm2] = estimate(cnm , snm , cnmcos , cnmsin , snmcos , snmsin , t0 ,t(2),n_cos );
[ ecnm3 ,esnm3] = estimate(cnm , snm , cnmcos , cnmsin , snmcos , snmsin , t0 ,t(3),n_cos );
[ ecnm4 ,esnm4] = estimate(cnm , snm , cnmcos , cnmsin , snmcos , snmsin , t0 ,t(4),n_cos );
 
h=450000 ; % the height for task 4
% compute the factors
for n=0:1:n_itsg
    anomaly(n+1) = GM/R^2*(n+1); % anomaly factor
    height (n+1) = (R/(R+h))^(n); % height factor
end

dist = GM/R ; %disturbance potential  factor
geoid = R ; % geoid hight factor

% ellipsoid factors
for n=0:1:n_grs
 anomalyg(n+1) = GM/R^2*(n+1);
 heightg (n+1) = (R/(R+h))^(n);
end

%cos , sin for all lambda for itsg
cosm = cosd([0:1:n_itsg]'*lambda);
sinm = sind([0:1:n_itsg]'*lambda);
%cos , sin for all lambda for ellipsiod
cosmg = cos([0:1:8]'*lambda);
sinmg = sin([0:1:8]'*lambda);
%fix the matrix dimension  for SHCs for ellipsiod
ccnmg= zeros(9,8);
ccnmg = horzcat(cnmg,ccnmg);
ssnmg= zeros(9,8);
ssnmg = horzcat(snmg,ssnmg);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:181 % loop over all thetas
    Pnm = Legendree(90-theta(i)',n_itsg); % all Legendre Functions for one theta itsg
    Pnmg = Legendree(90-theta(i)',n_grs);     % all Legendre Functions for one theta ellipsoid
    
    %Task 3
    % compute disturbance potential
    itsg_potetial (i,:)=  sum(dist*((cnm.*Pnm) * cosm + (snm.*Pnm) * sinm));% compute disturbance potential for itsg
    ellip_potential(i,:) =  sum(dist*((ccnmg.*Pnmg) * cosmg + (ssnmg.*Pnmg) * sinmg)); % compute disturbance potential for ellipsoid
    disturbance(i,:) = itsg_potetial(i,:) - ellip_potential(i,:); % substract the disturbance potential for ellip from itsg 
    
    % compute geoid height
    itsg_height (i,:)=  sum(geoid*((cnm.*Pnm) * cosm + (snm.*Pnm) * sinm));% compute geoid height for itsg
    ellip_height(i,:) =  sum(geoid*((ccnmg.*Pnmg) * cosmg + (ssnmg.*Pnmg) * sinmg)); % compute geoid height for ellipsoid
    geoid_height(i,:) = itsg_height(i,:) - ellip_height(i,:); % substract the geoid height for ellip from itsg
    
    % compute the gravity anomaly
    itsg_anomaly (i,:)=  anomaly*((cnm.*Pnm) * cosm + (snm.*Pnm) * sinm);% compute the gravity anomaly for itsg
    ellip_anomaly (i,:)=  anomalyg*((ccnmg.*Pnmg) * cosmg + (ssnmg.*Pnmg) * sinmg); % compute the gravity anomaly for ellipsoid
    gravity_anomaly(i,:) = (itsg_anomaly (i,:) - ellip_anomaly (i,:))*100000 ; % substract the gravity anomaly for ellip from itsg
    
    
    % compute disturbance potential at altitude 450 km
    itsg_potetial_h (i,:)=  height*(dist*((cnm.*Pnm) * cosm + (snm.*Pnm) * sinm)); %compute disturbance potential at altitude 450 km for itsg
    ellip_potential_h(i,:) =  heightg*(dist*((ccnmg.*Pnmg) * cosmg + (ssnmg.*Pnmg) * sinmg));  %compute disturbance potential at altitude 450 km for ellipsoid
    disturbance_h(i,:) = itsg_potetial_h(i,:) - ellip_potential_h(i,:); % substract the disturbance potential for ellip from itsg
    
    
    % compute geoid height for the new SHCs
    itsg_height_jan (i,:)=  sum(geoid*((ecnm1.*Pnm) * cosm + (esnm1.*Pnm) * sinm));
    itsg_height_april (i,:)=  sum(geoid*((ecnm2.*Pnm) * cosm + (esnm2.*Pnm) * sinm));
    itsg_height_july (i,:)=  sum(geoid*((ecnm3.*Pnm) * cosm + (esnm3.*Pnm) * sinm));
    itsg_height_october (i,:)=  sum(geoid*((ecnm4.*Pnm) * cosm + (esnm4.*Pnm) * sinm));
    
    % substract the geoid height for ellip from new SHCs
    geoid_height_jan(i,:) = itsg_height_jan(i,:) - ellip_height(i,:);
    geoid_height_april(i,:) = itsg_height_april(i,:) - ellip_height(i,:);
    geoid_height_july(i,:) = itsg_height_july(i,:) - ellip_height(i,:);
    geoid_height_october(i,:) = itsg_height_october(i,:) - ellip_height(i,:);

end

%comparing between the hights of the new SHCS and ITSG's SHCs  
diff_height1 = 1000*(geoid_height - geoid_height_jan);
diff_height2 = 1000*(geoid_height - geoid_height_april);
diff_height3 = 1000*(geoid_height - geoid_height_july);
diff_height4 = 1000*(geoid_height - geoid_height_october);


%
%Plotting Task 3-4
figure;

subplot(2,2,1);
 pcolor(L,T, disturbance);
 colorbar('fontsize',30), shading interp, daspect([1 1 1]);
 c= colorbar;
 c.Label.String = '(m^2/s^2)';
 title('Disturbance potential')
 caxis([-800 800])
hold on
plot(coast(:,1),(coast(:,2)),'g')
xlabel('Longitude'); ylabel('Latitude')
ax=gca;ax.XAxis.FontSize = 22;ax.YAxis.FontSize = 22;ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;ax.Title.FontSize = 24;ax.FontSize=18;
clear x

subplot(2,2,2);
 pcolor(L,T, gravity_anomaly);
 colorbar('fontsize',30), shading interp, daspect([1 1 1]);
 c= colorbar;
 c.Label.String = '(mGal)';
 title('Gravity anomalies')
 
hold on
plot(coast(:,1),(coast(:,2)),'g')
xlabel('Longitude'); ylabel('Latitude')
ax=gca;ax.XAxis.FontSize = 22;ax.YAxis.FontSize = 22;ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;ax.Title.FontSize = 24;ax.FontSize=18;clear x

subplot(2,2,3);
 pcolor(L,T, disturbance_h);
 colorbar('fontsize',30), shading interp, daspect([1 1 1]);
 c= colorbar;
 c.Label.String = ' (m^2/s^2)';
 title('Disturbance at altitude 450km')
 caxis([-800 800])
hold on
plot(coast(:,1),(coast(:,2)),'g')
xlabel('Longitude'); ylabel('Latitude')
ax=gca;ax.XAxis.FontSize = 22;ax.YAxis.FontSize = 22;ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;ax.Title.FontSize = 24;ax.FontSize=18;clear x

subplot(2,2,4);
 pcolor(L,T, geoid_height);
 colorbar('fontsize',30), shading interp, daspect([1 1 1]);
 c= colorbar;
 c.Label.String = '(m)';
 title('Geoid height')
caxis([-80 80]) 
hold on
plot(coast(:,1),(coast(:,2)),'g')
xlabel('Longitude'); ylabel('Latitude')
ax=gca;ax.XAxis.FontSize = 22;ax.YAxis.FontSize = 22;ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;ax.Title.FontSize = 24;ax.FontSize=18;clear x
% 

%
  
  



%%

%plotting the comparsion 
figure;

subplot(2,2,1);
 pcolor(L,T, diff_height1);
 colorbar('fontsize',30), shading interp, daspect([1 1 1]);
 c= colorbar;
 c.Label.String = '(mm)';
 title('Difference in geoid height between ITSG2018 (1st 2010) and in 1/1/2010')
 caxis([-5 5])
 hold on
plot(coast(:,1),(coast(:,2)),'g')
xlabel('Longitude'); ylabel('Latitude')
ax=gca;ax.XAxis.FontSize = 22;ax.YAxis.FontSize = 22;ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;ax.Title.FontSize = 24;ax.FontSize=18;
clear x
 
 subplot(2,2,2);
 pcolor(L,T, diff_height2);
 colorbar('fontsize',30), shading interp, daspect([1 1 1]);
 c= colorbar;
 c.Label.String = '(mm)';
 caxis([-5 5])
 title('Difference in geoid height between ITSG2018 (1st 2010) and in 1/4/2010')
 hold on
plot(coast(:,1),(coast(:,2)),'g')
xlabel('Longitude'); ylabel('Latitude')
ax=gca;ax.XAxis.FontSize = 22;ax.YAxis.FontSize = 22;ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;ax.Title.FontSize = 24;ax.FontSize=18;
clear x
 
 subplot(2,2,3);
 pcolor(L,T, diff_height3);
 colorbar('fontsize',30), shading interp, daspect([1 1 1]);
 c= colorbar;
 c.Label.String = '(mm)';
 caxis([-5 5])
 title('Difference in geoid height between ITSG2018 (1st 2010) and  in 1/7/2010')
 hold on
plot(coast(:,1),(coast(:,2)),'g')
xlabel('Longitude'); ylabel('Latitude')
ax=gca;ax.XAxis.FontSize = 22;ax.YAxis.FontSize = 22;ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;ax.Title.FontSize = 24;ax.FontSize=18;
clear x
 
 subplot(2,2,4);
 pcolor(L,T, diff_height4);
 colorbar('fontsize',30), shading interp, daspect([1 1 1]);
 c= colorbar;
 c.Label.String = '(mm)';
 title('Difference in geoid height between ITSG2018 (1st 2010) and in 1/10/2010')
caxis([-5 5])

hold on
plot(coast(:,1),(coast(:,2)),'g')
xlabel('Longitude'); ylabel('Latitude')
ax=gca;ax.XAxis.FontSize = 22;ax.YAxis.FontSize = 22;ax.XLabel.FontSize = 24;
    ax.YLabel.FontSize = 24;ax.Title.FontSize = 24;ax.FontSize=18;
clear x
