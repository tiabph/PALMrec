%% This script demonstrates the use of GPUgaussMLE
clear
clc
Nfits=1000;     %number of images to fit
bg=5;           %background fluorescence in photons/pixel/frame
Nphotons=300;   %expected photons/frame
Npixels=7;      %linear size of fit region in pixels. 
PSFsigma=1;     %PSF sigma in pixels

%generate a stack of images
coords=(Npixels-1)/2+zeros([Nfits 2]);
[out] = finiteGaussPSFerf(Npixels,PSFsigma,Nphotons,bg,coords);

%corrupt with Poisson noise 
data=noise(out,'poisson',1);
data=noise(data,'gaussian',1.5);
data=dip_array(data);

%fit and calculate speed
tic;
[X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(data,PSFsigma,20,2);
t=toc;
X=X-3;
Y=Y-3;
RMSE=sqrt(sum(X.^2+Y.^2)/1000);
fprintf('GPUgaussMLE has performed %g fits per second.\n',Nfits/t)

%report some details
s_x_found=std(X-coords(:,1));
meanCRLBx=mean(CRLBx);

fprintf('The standard deviation of x-position error is %g \n',s_x_found)
fprintf('The mean returned CRLB based x-position uncertainty is %g \n',meanCRLBx)


