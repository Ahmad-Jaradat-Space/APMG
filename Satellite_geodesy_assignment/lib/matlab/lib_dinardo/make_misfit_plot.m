function make_misfit_plot(swh,misfit,Mission_ID)

SWH=[0.05:0.1:8];

best_fit=zeros(size(SWH));

for j=1:length(SWH)
    
    
    best_fit(j)=nanmedian(misfit(swh>SWH(j)-0.1 & swh<SWH(j)+0.1));
    
    
end

figure('Position',get(0,'Screensize'));plot(swh,misfit,'o');hold on;h1=plot(SWH,best_fit,'r','LineWidth',3);
set(gca,'FontSize',15);title( [ Mission_ID ' 1Hz  SAR GOF' ]);grid on; xlabel([ 'SAR GPOD SWH [m]']);ylabel('Misfit between Model and Data [/]');axis tight;

[dum i0]=min(abs(SWH-2));
GOF=best_fit(i0)
