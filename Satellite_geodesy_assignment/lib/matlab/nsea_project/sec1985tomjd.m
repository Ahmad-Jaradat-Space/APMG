%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Diese Function convertiert ab einer bestimmten Referenzepoche die
%%%% die Daten nach MJD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=sec1985tomjd(sec)


%sec=338247747.1747291;
ref=date2mjd(1985,01,01,00,00,00);
refsec=ref*86400;
tage=fix(sec./86400);
secimtag=sec-(tage.*86400);
out=(ref+tage)+(secimtag./86400);

%mjd_sec=secimtag
%sekunden=refsec+sec;

