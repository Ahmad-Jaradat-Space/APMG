%%
%estimate SHCs

function [ ecnm ,esnm] = estimate(cnm , snm , cnmcos , cnmsin , snmcos , snmsin , t0 ,t,L )

estimated =[];
for n=0:1:L
    for m=0:1:n
        
ecnm(n+1,m+1)=[cnm(n+1,m+1)  ] + [cnmcos(n+1,m+1)]*cos(2*pi*((t-t0)/365.25)) + [cnmsin(n+1,m+1)]*sin(2*pi*((t-t0)/365.25));  
esnm(n+1,m+1)=[ snm(n+1,m+1) ] + [snmcos(n+1,m+1)]*cos(2*pi*((t-t0)/365.25)) + [snmsin(n+1,m+1)]*sin(2*pi*((t-t0)/365.25));
    end 
end
end