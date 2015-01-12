%% fit and plot
clear
clc
I=tiffread();
n=length(I(1,1,:));
[row col]=size(I(:,:,1));
fitl=(row-1)/2;
[Xd,Yd] = meshgrid(-fitl:fitl,-fitl:fitl);
xdata=[Xd,Yd];
V=[];
% ROIsum=zeros(row,col);
for i=1:n
    ROI=double(I(:,:,i));
    peak=double(ROI(fitl+1,fitl+1))-min(ROI(:));
    backg=min(ROI(:));
    maxm=max(ROI(:));
    [cy cx]=find(ROI==maxm);
    cy=round(mean(cy));
    cx=round(mean(cx));
    xp=[cx-fitl-1,cy-fitl-1,3,3,backg,peak,0];
    options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.001,'LargeScale','off');
    [lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian,xp,xdata,ROI,[],[],options);
    V(i,1)=lp(1);
    V(i,2)=lp(2);
    V(i,3)=lp(3);
    V(i,4)=lp(4);
    V(i,5)=lp(5);
    V(i,6)=lp(6);
%     ROI=ROI-lp(5);
%     Nm=sum(ROI(:))*10.4/5;
%     b=6*10.4/5;
%     a=160;
%     s=(lp(3)+lp(4))/2*160;
%     deltaR=sqrt((s^2+a^2/12)/Nm+8*pi*s^4*b^2/a^2/Nm^2);
%     xp=lp;
%     T=xp(7);
%     ROI=xp(6)*(exp(-0.5*((Xd-0)*cos(T)-(Yd-0)*sin(T)).^2./(xp(3)^2)-0.5*((Yd-0)*cos(T)+(Xd-0)*sin(T)).^2./(xp(4)^2)))+xp(5);
%     ROIsum=ROIsum+ROI;
    i
end
figure;
subplot(2,2,1);
plot(V(:,1),V(:,2),'.');
hold on
plot(V(1,1),V(1,2),'or')
plot(V(end,1),V(end,2),'og')
hold off
subplot(2,2,2)
D=[];
for i=1:n-1
    D(i)=sqrt((V(i+1,1)-V(1,1))^2+(V(i+1,2)-V(1,2))^2);
end
plot(D,'.');
hold on
plot(smooth(D,5),'r-')
axis([1 length(D) 0 2])
hold off
std(D(1:end))
std(D(1:end))*160
subplot(2,2,3)
hist(V(:,1),100);
subplot(2,2,4)
hist(V(:,2),100);


%% generate images
% peak=0.8;
% fitl=7;
% [Xd,Yd] = meshgrid(-fitl:fitl,-fitl:fitl);
% sx=2;
% I=zeros(15,15,10);
% for i=1:10
% cx=0.2*i-1;
% cy=0.2*i-1;
% ROI=peak*(exp(-0.5*(Xd-cx).^2./(sx^2)-0.5*(Yd-cy).^2./(sx^2)))+0.05;
% J=imnoise(ROI,'gaussian',0.1,0.2);
% J=double(J*1000);
% I(:,:,i)=J;
% end
% tiffwrite(I);

%% compute deltaR
% [row col]=size(I(:,:,1));
% xdata=[Xd,Yd];
% ROIsum=sum(I,3);
% ROI=ROIsum;
% peak=double(ROI(fitl+1,fitl+1))-min(ROI(:));
% backg=min(ROI(:));
% maxm=max(ROI(:));
% [cy cx]=find(ROI==maxm);
% cy=round(mean(cy));
% cx=round(mean(cx));
% xp=[cx-fitl-1,cy-fitl-1,2,2,backg,peak,0];
% options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',1000,'TolFun',0.001,'LargeScale','off');
% [lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian,xp,xdata,ROI,[],[],options);
% ROI=ROI-lp(5);
% Nm=sum(ROI(:))*10.4/5;
% b=6*10.4/5;
% a=160;
% s=(lp(3)+lp(4))/2*160
% deltaR=sqrt((s^2+a^2/12)/Nm+8*pi*s^4*b^2/a^2/Nm^2)
% xp=lp;
% % T=xp(7);
% % ROI=xp(6)*(exp(-0.5*((Xd-xp(2))*cos(T)-(Yd-xp(1))*sin(T)).^2./(xp(3)^2)-0.5*((Yd-xp(1))*cos(T)+(Xd-xp(2))*sin(T)).^2./(xp(4)^2)))+xp(5);
% ROI=xp(6)*(exp(-0.5*(Xd-xp(2)).^2./(xp(3)^2)-0.5*(Yd-xp(1)).^2./(xp(4)^2)))+xp(5);
% imshow(ROI,[])
