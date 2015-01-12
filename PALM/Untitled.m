% P=[];
% S=[];
% for i=1:1130
%     st=(i-1)*500+1;
%     et=i*500;
%     P(i)=mean(photonN(st:et,1));
%     S(i)=mean(sx(st:et,1));
% end
% plot(P,S)

clear
clc
peak=0.8;
fitl=4;
[Xd,Yd] = meshgrid(-fitl:fitl,-fitl:fitl);
xdata=[Xd,Yd];
sx=1.5;
cx=0;
cy=0;
ROI=peak*(exp(-0.5*(Xd-cx).^2./(sx^2)-0.5*(Yd-cy).^2./(sx^2)))+0.05;
I=zeros(9,9,1);
for i=1:1
J=imnoise(ROI,'gaussian',0.01,0.002);
J=imnoise(J,'poisson');
J=double(J*1000);
I(:,:,i)=J;
end
peak=max(I(:))-min(I(:));
backg=min(I(:));
xp=[0,0,1.5,1.5,backg,peak,0];
options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.01,'LargeScale','off');
[lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian,xp,xdata,I,[],[],options);
lp

%**************************************************************************
for i=1:9
    for j=1:9
       L(i,j)=(i-5)^2+(j-5)^2; 
    end
end
N=max(L(:));
for i=1:N
   II=(L==i);
   mm=mean(I(II));
   I(II)=mm;
end
peak=max(I(:))-min(I(:));
backg=min(I(:));
xp=[0,0,1.5,1.5,backg,peak,0];
options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.01,'LargeScale','off');
[lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian,xp,xdata,I,[],[],options);
lp
% tiffwrite(I);
% A=I(5,5,:);
% figure
% % plot(A(:))
% hold on
% A=I(4,4,:);
% plot(A(:),'b')
% A=I(6,6,:);
% plot(A(:),'r')
% A=I(4,6,:);
% plot(A(:),'g')
% A=I(6,4,:);
% plot(A(:),'y')