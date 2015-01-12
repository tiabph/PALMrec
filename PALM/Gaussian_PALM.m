function f=Gaussian_PALM(xp,xdata)
[a b]=size(xdata);
Xd=xdata(:,1:b/2);
Yd=xdata(:,b/2+1:b);
% T=xp(7);
f=xp(6)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(5);
% f=xp(5)*(exp(-0.5*((Xd-xp(2))*cos(T)-(Yd-xp(1))*sin(T)).^2./(xp(3)^2)-0.5*((Yd-xp(1))*cos(T)+(Xd-xp(2))*sin(T)).^2./(xp(4)^2)))+xp(6);