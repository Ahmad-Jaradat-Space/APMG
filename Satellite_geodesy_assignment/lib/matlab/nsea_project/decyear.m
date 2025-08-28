function output=decyear(ip)
% Reihenfolge Jahr / Monat / Tag / Stunde / Minute / Sekunde

output(:,1)=(ip(:,2)/12)+(ip(:,3)/365)+(ip(:,4)/8760)+(ip(:,5)/525600)+ip(:,1);
