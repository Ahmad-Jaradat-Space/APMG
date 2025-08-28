%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Umwandlung von MJD in Dezimales Jahr zwecks Diagrammdarstellung   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mjd=decyear2mjd(in)
%in=TUD_ALT(:,1);

out(:,1)=fix(in(:,1));                              % jahr
out(:,2)=in(:,1)-out(:,1);                          % Restjahr
out(:,3)=date2mjd(out(:,1),01,01,00,00,00);          % 1. Tag im Jahr
out(:,4)=date2mjd(out(:,1),12,31,00,00,00);          % Letzer Tag im Jahr
out(:,5)=out(:,4)-out(:,3)+1;                          % Tage im Jahr
out(:,6)=out(:,2).*out(:,5);                         % Tage bis aktuellem Datum 
out(:,7)=out(:,3)+out(:,6);                         % MJD bis zum aktuellen datum
%[out(:,8) out(:,9) out(:,10) out(:,11) out(:,12) out(:,13)]=mjd2date(out(:,7);
mjd(:,1)=out(:,7);
%mjd(:,2)=in(:,1);
%mjd(:,3)=mjd2decyear(out(:,7));


