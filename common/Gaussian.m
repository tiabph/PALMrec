function f=Gaussian(xp,xdata)
[a b]=size(xdata);
Xd=xdata(:,1:a);
Yd=xdata(:,a+1:2*a);
T=xp(7);
% f=xp(6)*(exp(-0.5*(Xd-xp(2)).^2./(xp(3)^2)-0.5*(Yd-xp(1)).^2./(xp(4)^2)))+xp(5);
f=xp(6)*(exp(-0.5*((Xd-xp(2))*cos(T)-(Yd-xp(1))*sin(T)).^2./(xp(3)^2)-0.5*((Yd-xp(1))*cos(T)+(Xd-xp(2))*sin(T)).^2./(xp(4)^2)))+xp(5);