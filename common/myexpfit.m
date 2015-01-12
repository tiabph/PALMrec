function [tau,x0,A,y0]=myexpfit(x,y,st)
% l=length(y);
y=y(:)';
xp=[5,st,(max(y)-min(y))*0.9,min(y)];
lb=[0.1,1,(max(y)-min(y))*0.2,0];
ub=[50,500,max(y)-min(y),max(y)];
options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.0001,'LargeScale','on');
[lp,resnorm,residual,exitflag]=lsqcurvefit(@oneexp,xp,x,y,lb,ub,options);
tau=lp(1);
x0=lp(2);
A=lp(3);
y0=lp(4);
% B=0;

function f=oneexp(xp,xdata)
f=xp(3)*(exp(-(xdata-xp(2))./xp(1)))+xp(4);

