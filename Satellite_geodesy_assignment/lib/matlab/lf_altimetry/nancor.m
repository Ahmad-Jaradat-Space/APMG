function cor=nancor(X,Y)

Co=nancov(X,Y);
cor=Co(1,2)/(sqrt(Co(1,1))*sqrt(Co(2,2)));