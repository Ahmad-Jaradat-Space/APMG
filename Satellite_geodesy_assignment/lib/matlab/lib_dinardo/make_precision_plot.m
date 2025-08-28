
function make_precision_plot(swh,std_SAR,std_RDS, type, unit,Mission_ID, RDS_LABEL)

SWH=[0.0:0.1:10];
best_fit=zeros(size(SWH));

for j=1:length(SWH)
    
    best_fit(j)=nanmedian(std_SAR(swh>SWH(j)-0.1 & swh<SWH(j)+0.1));
    
end



bbest_fit=zeros(size(SWH));

for j=1:length(SWH)
    
    bbest_fit(j)=nanmedian(std_RDS(swh>SWH(j)-0.1 & swh<SWH(j)+0.1));
    
end


figure('Position',get(0,'Screensize'));
%plot(swh,std_RDS,'bo');
hold on;plot(swh,std_SAR,'ro');
plot(SWH,best_fit,'k','LineWidth',3);
%plot(SWH,bbest_fit,'g','LineWidth',3);
set(gca,'FontSize',15);title(['1Hz ' type ' Precision']);grid on;
xlabel('SWH [m]');ylabel(['std ' type ' [' unit '] ']);axis tight;
%legend( [ Mission_ID ' ' RDS_LABEL ' POINTS' ], [Mission_ID ' SAR GPOD POINTS' ], [ Mission_ID ' SAR GPOD MEDIAN' ],[Mission_ID ' ' RDS_LABEL ' MEDIAN'])

[dum i0]=min(abs(SWH-2));
SAR_noise=best_fit(i0)
RDS_noise=bbest_fit(i0)
