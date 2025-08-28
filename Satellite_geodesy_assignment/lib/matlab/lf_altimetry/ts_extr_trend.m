function[ts_detrend,trend] = ts_extr_trend(date,ts)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE:        [ts_detrend,trend] = ts_extr_trend(date,ts)
%  
% DESCRIPTION:  fonction qui enleve le trend (signal lineaire) de chaque
%               serie temporelle  d une matrice de serie temporelles.
%
% INPUT:  
%   date:       vecteur des dates des series temporelles
%  
%   ts:         matrice de series temporelles en entree de format ligne k:
%               keme serie temporelle  
% 
% OUTPUT: 
%   ts_detrend: matrice des series temporelles en sortie de format ligne
%               k: keme serie temporelle sans trend
%
%  trend:       vecteur 2 colonnes avec le trend en premiere colonne et
%               l'ordonnee a l'origine dans la seconde pour chaque serie
%               temporelle   
%                        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  [lts,cts] = size(ts);

  trend = NaN*ones(lts,2);
  ts_detrend = NaN*ones(lts,cts);
  
  for i=1:lts
  
    %gestion des nan
    ts_tmp = ts(i,:);
    ind_nnan = not(isnan(ts_tmp));
    if (isempty(ind_nnan)==0)
      ts_tmp_nnan = ts_tmp(ind_nnan); 
      date_nnan = date(ind_nnan);
      
      %retrait de la tendance du signal
      
      [trend_tmp,cst,MU] = polyfit(date_nnan,ts_tmp_nnan,1);
      trend(i,1:2) = [trend_tmp(1)/MU(2),trend_tmp(2)-trend_tmp(1)*MU(1)/MU(2)];
      ts_detrend(i,:) = ts_tmp - (trend(i,1)*date + trend(i,2));
    else
      trend(i,1:2) = [NaN,NaN];
      ts_detrend (i,:) = ts_tmp;
    end
    
  end
    
