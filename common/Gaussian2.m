function f=Gaussian2(xp,xdata)
[a b]=size(xdata);
Xd=xdata(:,1:b/2);
Yd=xdata(:,b/2+1:b);
f=xp(1,5)*(exp(-0.5*(Xd-xp(1,1)).^2./(xp(1,3)^2)-0.5*(Yd-xp(1,2)).^2./(xp(1,4)^2)))+...
  xp(2,5)*(exp(-0.5*(Xd-xp(2,1)).^2./(xp(2,3)^2)-0.5*(Yd-xp(2,2)).^2./(xp(2,4)^2)))+xp(1,6);