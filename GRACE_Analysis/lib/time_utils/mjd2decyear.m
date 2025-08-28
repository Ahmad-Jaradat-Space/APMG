%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Umwandlung von MJD in Dezimales Jahr zwecks Diagrammdarstellung   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=mjd2decyear(in)
%in=date2mjd(2001,21,04,12,00,00);
%in=MTMW(:,1);

[yy mm dd h m s]=mjd2date(in);
in(:,2)=date2mjd(yy(:,1),01,01);
in(:,3)=date2mjd(yy(:,1),12,31);
in(:,4)=in(:,3)-in(:,2)+1;
in(:,5)=(in(:,1)-in(:,2));
in(:,6)=(in(:,5))./(in(:,4));
in(:,7)=yy(:,1)+in(:,6);
out=in(:,7);


%out=yy(:,1)+((mm(:,1)-1)/12)+((dd(:,1)-1)/365)+((h(:,1))/8760)+((m(:,1))/525600);




