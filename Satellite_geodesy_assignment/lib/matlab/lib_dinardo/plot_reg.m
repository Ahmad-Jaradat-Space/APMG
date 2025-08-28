function [brob]=plot_reg(x,y)

brob = robustfit(x,y);
%x=[ones(size(x)); x];
%[brob] = regress(y.',x.');
yf=brob(1)+brob(2).*x;

yf(yf<min(y))=nan;
yf(yf>max(y))=nan;
plot(x,yf,'r','LineWidth',2,'Parent',gca);

end

