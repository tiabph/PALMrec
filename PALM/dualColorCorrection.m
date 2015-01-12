clear
clc
pathname = 'G:\lzl\20140104\4+';
str = strcat(pathname,'\488.tif');
I1 = tiffread(str);
fitl = (size(I1,1)-1)/2;
[X Y]= meshgrid(-fitl:fitl,-fitl:fitl);
xdata=[X Y];
Fd = double(I1);
peak = max(Fd(:));
bg = min(Fd(:));
xp=[0,0,1.5,1.5,bg,peak];
options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.01,'LargeScale','on');
[lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian,xp,xdata,Fd,[],[],options);
V1(1,1) = lp(2);
V1(1,2) = lp(1);

str = strcat(pathname,'\561.tif');
I2 = tiffread(str);
Fd = double(I2);
peak = max(Fd(:));
bg = min(Fd(:));
xp=[0,0,1.5,1.5,bg,peak];
options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.01,'LargeScale','on');
[lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian,xp,xdata,Fd,[],[],options);
V2(1,1) = lp(2);
V2(1,2) = lp(1);

V = V2-V1;