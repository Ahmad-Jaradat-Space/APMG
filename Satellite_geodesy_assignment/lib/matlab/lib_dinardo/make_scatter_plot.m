
function make_scatter_plot(x,y,type1,type2,unit,No)

nump2=sum(not(isnan(x-y)));
perc2=sum(not(isnan(x-y)))*100/No;


figure('Position',get(0,'Screensize'));

if nump2==0
    
    disp('No point to make scatterplot--> scatter plot skipped')
    return
end

scattercloud(x,y,50,0.75,'none',colormap(flipud(hot(64))));grid on;
set(gca,'FontSize',20);
xlabel([ type1 ' [' unit ']'],'Interpreter','none');ylabel([ type2 ' [' unit ']'],'Interpreter','none');
hold on;plot(x,x,'k','LineWidth',3);

scatt_stat=plot_reg(x,y);
set(findobj(gcf,'type','axes'),'layer','top');
set(findobj(gcf,'type','axes'),'Xtickmode','auto');
set(findobj(gcf,'type','axes'),'Ytickmode','auto')
scatt_stat;colorbar;
text(0.5,0.35,['Regression Slope: ' num2str(scatt_stat(2),'% 10.3f')],'units','normalized','FontSize',19);
text(0.5,0.30,['STDD: ' num2str(nanstd(x-y),'% 10.3f') unit],'units','normalized','FontSize',19); 
text(0.5,0.25,['Bias: ' num2str(nanmedian(y-x),'% 10.3f') unit],'units','normalized','FontSize',19);
text(0.5,0.20,['NP : ' num2str(nump2,'% 10.3d')  ],'units','normalized','FontSize',19);
%text(0.5,0.12,['NP retained(%): ' num2str(perc2,'% 10.3f'),' % of ',  num2str(No,'% 10.3d')  ],'units','normalized','FontSize',19);
text(1.16,0.999,'%','units','normalized','FontSize',22);set(gca,'FontSize',25);set(gca,'FontSize',20);