function fit
str='F:\drift\bead1.tif';
I=tiffread(str,[1 5000]);
pause(eps)
n=length(I);
wid=size(I,1);
[Xd,Yd] = meshgrid(1:wid,1:wid);
xdata=[Xd,Yd];
B=zeros(n,2);
for i=1:n
    Fd=double(I(:,:,i));
    backg=min(Fd(:));
    peak=max(Fd(:));
    [cy cx]=find(Fd==peak);
    cy=mean(cy);
    cx=mean(cx);
    xp=[cx,cy,2,2,backg,peak];
    options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',1000,'TolFun',0.01,'Algorithm','levenberg-marquardt');
    [lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian_PALM1,xp,xdata,Fd,[],[],options);
    B(i,2)=lp(1);
    B(i,1)=lp(2);
%     [yy xx]=weight_centrid(Fd,wid);
%     if abs(yy-cy)>1 || abs(xx-cx)>1
%         yy=cy;
%         xx=cx;
%     end
%     if abs(B(i,1)-yy)>2 || abs(B(i,2)-xx)>2
%         B(i,1)=yy;
%         B(i,2)=xx;
%     end
    i
end

% for i=2:length(B)-1
%     if abs(B(i,1)-B(i-1,1))>3 && abs(B(i-1,1)-B(i+1,1))<1
%         B(i,1)=(B(i-1,1)+B(i+1,1))/2;
%     end
%     if abs(B(i,2)-B(i-1,2))>3 && abs(B(i-1,2)-B(i+1,2))<1
%         B(i,2)=(B(i-1,2)+B(i+1,2))/2;
%     end
% end

% str='F:\mEos3\2011-9-29\WT\beads.tif';
% I=tiffread(str);
% Fd=double(I);
% backg=min(Fd(:));
% peak=max(Fd(:));
% [cy cx]=find(Fd==peak);
% cy=mean(cy);
% cx=mean(cx);
% xp=[cx,cy,2,2,backg,peak];
% options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',1000,'TolFun',0.01,'Algorithm','levenberg-marquardt');
% [lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian_PALM1,xp,xdata,Fd,[],[],options);
% x=lp(1);
% y=lp(2);
% [yy xx]=weight_centrid(Fd,wid);
% if abs(yy-cy)>1 || abs(xx-cx)>1
%     yy=cy;
%     xx=cx;
% end
% if abs(y-yy)>2 || abs(x-xx)>2
%     y=yy;
%     x=xx;
% end
% B(:,1)=B(:,1)-y;
% B(:,2)=B(:,2)-x;
figure;plot(B(:,1),'b');hold on;plot(B(:,2),'r');

function f=Gaussian_PALM1(xp,xdata)
[a b]=size(xdata);
Xd=xdata(:,1:b/2);
Yd=xdata(:,b/2+1:b);
f=xp(6)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(5);

function [y,x]=weight_centrid(ROI,w)
sumx=0;
sumy=0;
for i=1:w
    for j=1:w
        sumx=sumx+ROI(i,j)*j;
        sumy=sumy+ROI(i,j)*i;
    end
end
y=sumy/sum(ROI(:));
x=sumx/sum(ROI(:));