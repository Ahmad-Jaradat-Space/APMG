function [time_out,lon_out,lat_out,ts_out] = lec_ts_rep_netcdf_duacs(rep_in,ext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE:       [time_out,lon_out,lat_out,ts_out] = lec_ts_rep_netcdf_duacs(rep_in)
%  
% DESCRIPTION: Fonction qui lit les fichiers contenant des cartes de
%              series temporelles de ssh en netcdf d aviso delivree par
%              cls a l adresse suivante:
%              ftp://ftpsedr.cls.fr/pub/oceano/AVISO/SSH/duacs/regional-mfstep/dt/ref/msla/merged/h/
%
% INPUT:
%   rep_in:    repertoire contenant tous les fichiers netcdf de ssh
%              delivres par AVISO. rep_in doit se teminer par /.
%   ext:       extension des fichiers que l on veut lire: * pour tous les
%              fichiers ou *qd* pour tous les fichiers qui contiennent qd
%              dans leur nom.   
%
% OUTPUT:  
%   lon_out:   Vecteur colonne des longitudes. les longitudes sont modulo
%              360 deg.
%   lat_out:   Vecteur colonne des latitudes. 
%   ts_out:    (mm) Matrice dont les lignes sont les series temporelles de
%              ssh. Format:  lig k : serie temporelle de ssh du point de
%              longitude lon_out(k) et de latitude lat_out(k).
%              Les valeurs abberrentes sont supposees de valeurs > 1000
%              dans les fichiers input et sont egalee a NaN dans ts_out.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %recuperation des noms des fichiers dans rep_in
  %scale=0.0001% scaling factor to transform in m
  scale=0.1% scaling factor to transform in mm
  %fillvalue = -2147483647 ;
  fillvalue = -2147483647 ;
  rep_st = dir([rep_in,ext]);
  [nrep_st, dum] = size(rep_st);

  pf = 1;
  nb_fic = nrep_st;
  if (ext == '*') %(on enleve rep_in/. et rep_in/.. quand ils sont comptes)
    pf = 3;
    nb_fic = nrep_st-2;
  end
  
  %cas du premier fichier
  %nom du fichier
  fic_in = strcat(rep_in,rep_st(pf).name);
    
  %ouverture du fichier netcdf
  ncid = netcdf.open(fic_in, 'nowrite');  

  %recuperation des var id
  id_time = netcdf.inqVarID(ncid,'time');
  id_lat = netcdf.inqVarID(ncid,'lat');
  id_lon = netcdf.inqVarID(ncid,'lon');
  id_grid = netcdf.inqVarID(ncid,'sla');
  
  %recuperation des valeurs
  time = netcdf.getVar(ncid,id_time);
  lat = netcdf.getVar(ncid,id_lat);
  lon = netcdf.getVar(ncid,id_lon);
  %grid = transpose(netcdf.getVar(ncid,id_grid));
  grid = netcdf.getVar(ncid,id_grid);
  
  %ecriture des variables de sorties
  [llat,clat] = size(lat);  
  [llon,clon] = size(lon);  
  lon_out = repmat(lon,llat,1);
  lat_out = repmat(lon,llat,1);
  lat_out(1:llon,1) = lat(1);
  % allocate matrix and vector
  ts_out = repmat(NaN,llat*llon,nb_fic);
  time_out = repmat(NaN,1,nb_fic);
  ts_out(1:llon,1) = (grid(:,1)); %les donnees DUACS sont en 0.1 mm
  %time_out(1) =  time; %les donnees CLS sont en cm
  jd1950 = date2jd(1950);
  jd1985  = date2jd(1985);
  tsec = (time - jd1985 + jd1950)*86400;
%   tmjd = sec1985_to_mjd(tsec)
%   yfrac=mjd2decyear(tmjd);
  time_out(1) = tsec;
  %    
  for k=2:llat
  lat_out(((k-1)*llon+1):(k*llon),1) = lat(k);
  ts_out(((k-1)*llon+1):(k*llon),1) = (grid(:,k)); %les donnees DUACS sont en 0.1 mm
  end
  
  %on ferme le fichier netcdf
  netcdf.close(ncid);
  
  %on boucle sur les fichiers du repertoire (on enleve bien sur
  %rep_in/. et rep_in/.. quand ils sont comptes)
  if (nrep_st > pf)
  
    for i = (pf+1):nrep_st
  
      compt = i;
      if (pf==3)
        compt = i-2;
      end
      
      %nom du fichier
      fic_in = strcat(rep_in,rep_st(i).name);
      
      %ouverture du fichier netcdf
      ncid = netcdf.open(fic_in, 'nowrite');  

      %recuperation des valeurs
      time_current = netcdf.getVar(ncid,id_time);
      lat_current = netcdf.getVar(ncid,id_lat);
      lon_current = netcdf.getVar(ncid,id_lon);
      %grid = transpose(netcdf.getVar(ncid,id_grid));
      grid = netcdf.getVar(ncid,id_grid);
            
      %ecriture des variables de sorties
            
      [llat_current,clat_current] = size(lat_current);  
      [llon_current,clon_current] = size(lon_current);
      if ( ([llat_current,clat_current] ~= [llat,clat]) | ([llon_current,clon_current] ~= [llon,clon]) )
        error(strcat('il n y a pas le meme nombre de longitudes ou latitude dans les fichiers netcdf du repertoire ',rep_in));
      end
      if ( any(lat_current-lat) | any(lon_current-lon) )
        error(strcat('il n y a pas les memes valeurs de longitude ou latitude dans les fichiers netcdf du repertoire ',rep_in));
      end
      ts_out(1:llon,compt) =  (grid(:,1)); %les donnees CLS sont en cm
      % convert time in MODEL from days1950 to sec1985,mjd and yfrac
      time=time_current;
      jd1950 = date2jd(1950);
      jd1985  = date2jd(1985);
      tsec = (time - jd1985 + jd1950)*86400;
      time_out(compt)=tsec;
%       tmjd = sec1985_to_mjd(tsec);
%       yfrac=mjd2decyear(tmjd);
%       time_out(compt) = yfrac;
      
      for k=2:llat
        ts_out(((k-1)*llon+1):(k*llon),compt) =  (grid(:,k)); %les donnees CLS sont en cm
      end

      %on ferme le fichier netcdf
      netcdf.close(ncid);


    end      
  end
  

  %on met les valeurs abberentes a NaN et les longitude a modulo 360
  %f1 = find(ts_out > 1000);
  f1 = find(ts_out < fillvalue+1);
  ts_out(f1) = NaN;
  ts_out =ts_out*scale; %transform in mm 
  f2 = find(lon_out > 360.| lon_out < 0.);
  lon_out(f2) = mod(lon_out(f2),360.);
 
