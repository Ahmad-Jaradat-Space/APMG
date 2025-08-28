function [time_out,ts_out]=moy_mens2an(time_in,ts_in)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Fonction qui moyenne les valeurs mensuelles en valeurs annuelles.
% La moyennation est faite annee par annee. On consid?re que si
% il y a k ou plus mois par annee cela est suffisant pour calculer
% la moyenne annuelle.
% input: ts_in est la matrice des series temporelles mensuelles en entree
%        time_in le vecteur temporel correspondant en annees.ATTENTION: Il doit
%        avoir un nombre entier d annee et les mois doivent
%        etre datees du milieu du mois: par exemple le 1 ere
%        mois de 1992 correspond a la date 1992+0.5/12 = 1.992041666666667e+03
%  
% changed due to error atline 37 column 49 (commented)
%
% output: ts_out est la matrice des series temporelles annuelles en sortie
%        time_out le vecteur temporel correspondant en annees  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [lts_in,cts_in] = size(ts_in);
  
  % seuil de difference entre deux doubles pour les fonctions
  % intersect_dbl et union_dbl ismember_dbl
  epsilon = 1.e-8;
  
 % recuperation du nombre d annees et verification qu il y a un
 % nombre entier d annee dans les donnees et que le vecteur date
 % time_in commence bien au debut d une annee
 nb_an = fix(cts_in/12);
 test_nb_an = rem(cts_in,12);
 time_in_min = min(time_in);
 test_time_in_min = (time_in_min - 0.5/12) - floor(time_in_min); 
 
 if ( (test_nb_an ~= 0) )
   error('il n y a pas un nombre entier dannee dans les donnees a moyenner');
 end
 
%  if (not(ismember(test_time_in_min,[0],epsilon))
%    error('les donnees a moyenner ne commencent pas au debut dune annee ');
%  end
 
 
 % calcul de la moyenne annuelle
 %pour chaque annee qui contient k ou plus de semaines on calcule la moyenne
 k=9;
 
 
 time_out = NaN(1,nb_an);
 ts_out = NaN(lts_in,nb_an);
 
 time_tmp = NaN(1,12);
 ts_tmp = NaN(lts_in,12);
 denom = NaN(lts_in,1);
 
 for i = 1:nb_an
 
   time_tmp = time_in(1+12*(i-1):12*i);
   time_out(i) = floor(time_tmp(1)) + 0.5;
   
   ts_tmp = ts_in(:,1+12*(i-1):12*i);
   f=find(isnan(ts_tmp)==1);
   ind = ones(size(ts_tmp));
   ts_tmp(f) = 0;
   ind(f) = 0;
   denom = sum(ind,2);   
   %pour chaque annee qui contient k ou plus mois on calcule la moyenne
   %annuelle
   f1 = find(denom>=k);
   ts_out(f1,i) = sum(ts_tmp(f1,:),2)./denom(f1);
   
 end
 
