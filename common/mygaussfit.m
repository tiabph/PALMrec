function [sigma,mu,A,B]=mygaussfit(x,y,minStd,maxStd)
l=length(y);
y=y(:)';%-min(y);
xp=[10,(l+1)/2,max(y)-min(y),min(y)];
lb=[minStd,(l+1)/2-3,0,0];
ub=[maxStd,(l+1)/2+3,max(y),max(y)];
options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',1000,'TolFun',0.01,'LargeScale','on');
[lp,resnorm,residual,exitflag]=lsqcurvefit(@onegauss,xp,x,y,lb,ub,options);
sigma=lp(1);
mu=lp(2);
A=lp(3);
B=lp(4);
% B=0;

function f=onegauss(xp,xdata)
f=xp(3)*(exp(-0.5*(xdata-xp(2)).^2./(xp(1)^2)))+xp(4);

