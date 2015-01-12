clear
clc
for k=5:5
bg=5;           % background photon     
N=100*k;        % total photon
sigma=2;      % PSF sigma
% sigma=1.3+0.2*randn(1000,1);
fitl=11;         % fit size
n=1000;         % image number
[Xd,Yd] = meshgrid(1:fitl,1:fitl);
sx=sigma;
sy=sigma;
cx=(fitl+1)/2;
cy=(fitl+1)/2;
ROI=(exp(-0.5*(Xd-cx).^2./(sx^2)-0.5*(Yd-cy).^2./(sy^2)));
peak=N/sum(ROI(:));
ratio=peak/0.5;
b=bg/peak*0.5;
I=zeros(fitl,fitl,n);
ROI=0.5*(exp(-0.5*(Xd-cx).^2./(sx^2)-0.5*(Yd-cy).^2./(sy^2)))+b;    
for i=1:n
% sx=sigma(i);
% sy=sigma(i);
% ROI=(exp(-0.5*(Xd-cx).^2./(sx^2)-0.5*(Yd-cy).^2./(sy^2)));
% peak=N/sum(ROI(:));
% ratio=peak/0.5;
% b=bg/peak*0.5;
% ROI=0.5*(exp(-0.5*(Xd-cx).^2./(sx^2)-0.5*(Yd-cy).^2./(sy^2)))+b;    
J=imnoise(ROI,'gaussian',0,0.001);
J=uint16(J*ratio);
J=imnoise(J,'poisson');
I(:,:,i)=J;
end
% tiffwrite(uint16(I),['N',num2str(N),'.tif']);
end
data=I;

% clear
% clc
% data=tiffread();
n=length(data(:,:,1));
tic
[X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(single(data),sigma,10,2);
X=X+1;
toc
RMSE_position=sqrt(sum((X-4).^2)/n)
MAE_position=mean(abs(X-4))
RMSE_sigma=sqrt(sum((S-1).^2)/n)
MAE_sigma=mean(abs(S-1))
RMSE_N=sqrt(sum((N-500).^2)/n)
MAE_N=mean(abs(N-500))
RMSE_BG=sqrt(sum((BG-5).^2)/n)
MAE_BG=mean(abs(BG-5))
