
function make_npoint_distance_plot_2(DISTANCE,RDS,SAR,MOD,Mission_ID,Regional_Model_ID,RDS_LABEL)

dist=[0:200:10e3];
xx = dist(1:end-1) + 100;
y1 = nan(1,length(dist)-1);
y2 = nan(1,length(dist)-1);
y3 = nan(1,length(dist)-1);
for i=1:length(dist)-1
    y1(i) = nnz(~isnan(SAR( DISTANCE > dist(i) & DISTANCE < dist(i+1))));
    y2(i) = nnz(~isnan(RDS(DISTANCE > dist(i) & DISTANCE < dist(i+1))));
    y3(i) = nnz(~isnan(MOD(DISTANCE > dist(i) & DISTANCE < dist(i+1))));
end
numMODi=sum(y3);
numSLApSARi=sum(y1);
numSLApRDSi=sum(y2);

figure('Position',get(0,'Screensize'))
bar(xx,y2,'histc');
hold on
bar(xx,y1,'histc');
bar(xx,y3,'histc');
hl=legend([Mission_ID ' ' RDS_LABEL ] , [Mission_ID ' SAR' ] ,[Regional_Model_ID ' OCEAN MODEL']);
set(hl,'Interpreter','none');
h = findobj(gca,'Type','patch');
set(gcf,'Renderer','opengl')
set(h(1),'EdgeColor','k');hold on;
set(h(1),'FaceColor','g');set(gca,'FontSize',20)
set(h(2),'FaceColor','r');set(gca,'FontSize',20)
set(h(3),'FaceColor','b');set(gca,'FontSize',20)
set(gca,'FontSize',15);title('Number of 20Hz points in 200 meter step vs DISTANCE to COAST','FontSize',20);grid on;
xlabel('DISTANCE to coast [m]','FontSize',20);
ylabel('NPoints 20Hz','FontSize',20);axis tight
text(0.8,0.80,['NP SAR GPOD: ' num2str(numSLApSARi,'% 10.3d')],'units','normalized','FontSize',20);
text(0.8,0.75,['NP ' RDS_LABEL ' : ' num2str(numSLApRDSi,'% 10.3d')],'units','normalized','FontSize',20);
text(0.8,0.70,['NP MODEL: ' num2str(numMODi,'% 10.3d')],'units','normalized','FontSize',20);
