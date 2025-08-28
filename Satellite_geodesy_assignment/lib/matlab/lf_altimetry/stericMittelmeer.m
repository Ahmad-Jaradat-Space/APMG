

dataT = nc_reader('/home/fenoglio/softw_ex/IGG_thermoSter_THETA_MONMEAN.nc');
dataS = nc_reader('/home/fenoglio/softw_ex/IGG_haloSter_SALT_MONMEAN.nc');

%%
mask = ~isnan(dataT.thermoSter(:,:,1,1)); % & data.lon>-6;

id = isnan(dataT.thermoSter);
temp = dataT.thermoSter;
temp(id)=0;
salt = dataS.haloSter;
salt(id)=0;

ster=[];
mm2 = [];
for j=1:12
    
    tsm = temp(:,:,:,j);
    tsm = sum(tsm,3);
   
    ssm = salt(:,:,:,j);
    ssm = sum(ssm,3);
    hv = tsm+ssm;
    ster = [ ster hv(:)]; % neue spalte
    
    hv = hv(:);
    mm2 = [mm2 nanmean(hv(mask(:)))];
    
end


%%
mster = nanmean(ster(:,119:1919),2); % Mittelwert in jedem Punkt ueber alle Monate
ster = ster -repmat(mster,1,size(ster,2)); % sterische anomalie in jedem punkt

mm =[];

for j=1:size(ster,2)
    mm=[mm nanmean(ster(mask(:),j))];

end

figure, plot(mm); hold on; plot(mm2+(mm(1)-mm2(1)),'r');
% blau (mm) ist das was man will
