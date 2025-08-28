function RMSval=rms(S)

[m n]=size(S);
for i=1:n
f=find(isnan(S(:,i))==0);

RMSval(i)=sqrt(sum(S(f,i).^2)/length(S(f,i)));
end