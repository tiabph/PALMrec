function VH=PALM_fitting_iPALM(VH)

gain=VH.parameter.fitting.gain;
factor=VH.parameter.fitting.factor;
sigma=VH.parameter.fitting.sigma;
a=VH.parameter.fitting.pixelsize;
fitsize=VH.parameter.fitting.fitsize;
fitl=(VH.parameter.fitting.fitsize-1)/2;
photonflag=VH.parameter.fitting.photonflag;
minSigma=VH.parameter.reconstruction.minSigma;
maxSigma=VH.parameter.reconstruction.maxSigma;
gap=VH.parameter.linking.gap;

row=VH.row;
column=VH.column;
m=length(VH.DV);
% fitl=round(3*sigma);

set(VH.text10,'string','Generating ROI ...')
mywaitbar(0,VH.axes2,'');
ROIfit=single(zeros(2*fitl+1,2*fitl+1,m));
[Xd Yd]= meshgrid(-fitl:fitl,-fitl:fitl);
for i=1:m
    n=length(VH.DV(i,1).trackInfo(:,1));
    cx=round(VH.DV(i,1).trackInfo(1,3));
    cy=round(VH.DV(i,1).trackInfo(1,2));
    up=cx-fitl; 
    bottom=cx+fitl; 
    left=cy-fitl;
    right=cy+fitl;   
    if up>=1 && left>=1 && bottom<=row  && right<=column  && n<VH.ImageNumber             
        frame=VH.DV(i,1).trackInfo(:,1); 
        I=double(VH.A(up:bottom,left:right,frame));
        ROI=sum(I,3);
        ROI=ROI-min(ROI(:));  
        ROIfit(:,:,i)=single(ROI);            
    end 
    if (floor(i/m*100))>(floor((i-1)/m*100))
    mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
end

fitInfo=[];
set(VH.text10,'string','GPU fitting ...')
N=floor(m/10000);
for i=1:N+1
    st=(i-1)*10000+1;
    et=i*10000;
    et=min(et,m);
    ROI=ROIfit(:,:,st:et)*factor/gain;
    tic
    [X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(single(ROI),sigma,20,2);
    toc
    fitInfo(st:et,1)=X+1;
    fitInfo(st:et,2)=Y+1;
    fitInfo(st:et,3)=N;
    fitInfo(st:et,4)=BG;
    fitInfo(st:et,5)=S;
end

L=length(fitInfo(:,5));
CX=zeros(L,4);
CY=zeros(L,4);
%**************************************************************************
filebase=VH.filebase(1:end-4);
str=strcat(VH.pathname,filebase,'_405.tif');
A1=tiffread(str,[1 VH.ImageNumber]);
set(VH.text10,'string','Generating ROI ...')
mywaitbar(0,VH.axes2,'');
ROIfit1=single(zeros(2*fitl+1,2*fitl+1,m));
for i=1:m
    n=length(VH.DV(i,1).trackInfo(:,1));
    cx=round(VH.DV(i,1).trackInfo(1,3));
    cy=round(VH.DV(i,1).trackInfo(1,2));
    up=cx-fitl; 
    bottom=cx+fitl; 
    left=cy-fitl;
    right=cy+fitl;   
    if up>=1 && left>=1 && bottom<=row  && right<=column  && n<VH.ImageNumber             
        frame=VH.DV(i,1).trackInfo(:,1); 
        I=double(A1(up:bottom,left:right,frame));
        ROI=sum(I,3);
%         ROI=ROI-min(ROI(:));  
        ROIfit1(:,:,i)=single(ROI);            
    end 
    if (floor(i/m*100))>(floor((i-1)/m*100))
    mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
end

fitInfo1=[];
set(VH.text10,'string','GPU fitting ...')
N=floor(m/10000);
for i=1:N+1
    st=(i-1)*10000+1;
    et=i*10000;
    et=min(et,m);
    ROI=ROIfit1(:,:,st:et);%*factor/gain;
    tic
    [X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(single(ROI),sigma,20,2);
    toc
%     num = length(X);
%     M3 = zeros(num,1);
%     for j=1:num
%         ccx = X(j,1)-fitl;
%         ccy = Y(j,1)-fitl;
%         sx = median(S);
%         ROI1 = (exp(-0.5*(Xd-ccx).^2./(sx^2)-0.5*(Yd-ccy).^2./(sx^2)));
%         Dist = ((Xd-ccx).^2+(Yd-ccy).^2).^(3/2);
%         ROI1 = ROI1.*Dist;
%         ROI1 = ROI1/sum(ROI1(:));
%         ROI1 = ROI1.*single(ROI(:,:,j));
%         ROI1 = ROI1.*(ROI1>0);
%         M3(j,1) = sum(ROI1(:));
%     end
%     fitInfo2(st:et,2)=M3;
    fitInfo1(st:et,1)=N;
    CX(st:et,1)=X+1;
    CY(st:et,1)=Y+1;
end
%**************************************************************************

%**************************************************************************
filebase=VH.filebase(1:end-4);
str=strcat(VH.pathname,filebase,'_473.tif');
A1=tiffread(str,[1 VH.ImageNumber]);
set(VH.text10,'string','Generating ROI ...')
mywaitbar(0,VH.axes2,'');
ROIfit1=single(zeros(2*fitl+1,2*fitl+1,m));
for i=1:m
    n=length(VH.DV(i,1).trackInfo(:,1));
    cx=round(VH.DV(i,1).trackInfo(1,3));
    cy=round(VH.DV(i,1).trackInfo(1,2));
    up=cx-fitl; 
    bottom=cx+fitl; 
    left=cy-fitl;
    right=cy+fitl;   
    if up>=1 && left>=1 && bottom<=row  && right<=column  && n<VH.ImageNumber             
        frame=VH.DV(i,1).trackInfo(:,1); 
        I=double(A1(up:bottom,left:right,frame));
        ROI=sum(I,3);
%         ROI=ROI-min(ROI(:));  
        ROIfit1(:,:,i)=single(ROI);            
    end 
    if (floor(i/m*100))>(floor((i-1)/m*100))
    mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
end

set(VH.text10,'string','GPU fitting ...')
N=floor(m/10000);
for i=1:N+1
    st=(i-1)*10000+1;
    et=i*10000;
    et=min(et,m);
    ROI=ROIfit1(:,:,st:et);%*factor/gain;
    tic
    [X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(single(ROI),sigma,20,2);
    toc
%     num = length(X);
%     M3 = zeros(num,1);
%     for j=1:num
%         ccx = X(j,1)-fitl;
%         ccy = Y(j,1)-fitl;
%         sx =  median(S);
%         ROI1 = (exp(-0.5*(Xd-ccx).^2./(sx^2)-0.5*(Yd-ccy).^2./(sx^2)));
%         Dist = ((Xd-ccx).^2+(Yd-ccy).^2).^(3/2);
%         ROI1 = ROI1.*Dist;
%         ROI1 = ROI1/sum(ROI1(:));
%         ROI1 = ROI1.*single(ROI(:,:,j));
%         ROI1 = ROI1.*(ROI1>0);
%         M3(j,1) = sum(ROI1(:));
%     end
%     fitInfo2(st:et,1)=M3;
    fitInfo1(st:et,2)=N;
    CX(st:et,2)=X+1;
    CY(st:et,2)=Y+1;
end
%**************************************************************************

%**************************************************************************
filebase=VH.filebase(1:end-4);
str=strcat(VH.pathname,filebase,'_561.tif');
A1=tiffread(str,[1 VH.ImageNumber]);
set(VH.text10,'string','Generating ROI ...')
mywaitbar(0,VH.axes2,'');
ROIfit1=single(zeros(2*fitl+1,2*fitl+1,m));
for i=1:m
    n=length(VH.DV(i,1).trackInfo(:,1));
    cx=round(VH.DV(i,1).trackInfo(1,3));
    cy=round(VH.DV(i,1).trackInfo(1,2));
    up=cx-fitl; 
    bottom=cx+fitl; 
    left=cy-fitl;
    right=cy+fitl;   
    if up>=1 && left>=1 && bottom<=row  && right<=column  && n<VH.ImageNumber             
        frame=VH.DV(i,1).trackInfo(:,1); 
        I=double(A1(up:bottom,left:right,frame));
        ROI=sum(I,3);
%         ROI=ROI-min(ROI(:));  
        ROIfit1(:,:,i)=single(ROI);            
    end 
    if (floor(i/m*100))>(floor((i-1)/m*100))
    mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
end

set(VH.text10,'string','GPU fitting ...')
N=floor(m/10000);
for i=1:N+1
    st=(i-1)*10000+1;
    et=i*10000;
    et=min(et,m);
    ROI=ROIfit1(:,:,st:et);%*factor/gain;
    tic
    [X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(single(ROI),sigma,20,2);
    toc
%     num = length(X);
%     M3 = zeros(num,1);
%     for j=1:num
%         ccx = X(j,1)-fitl;
%         ccy = Y(j,1)-fitl;
%         sx =  median(S);
%         ROI1 = (exp(-0.5*(Xd-ccx).^2./(sx^2)-0.5*(Yd-ccy).^2./(sx^2)));
%         Dist = ((Xd-ccx).^2+(Yd-ccy).^2).^(3/2);
%         ROI1 = ROI1.*Dist;
%         ROI1 = ROI1/sum(ROI1(:));
%         ROI1 = ROI1.*single(ROI(:,:,j));
%         ROI1 = ROI1.*(ROI1>0);
%         M3(j,1) = sum(ROI1(:));
%     end
%     fitInfo2(st:et,4)=M3;
    fitInfo1(st:et,3)=N;
    CX(st:et,3)=X+1;
    CY(st:et,3)=Y+1;
end
%**************************************************************************

%**************************************************************************
filebase=VH.filebase(1:end-4);
str=strcat(VH.pathname,filebase,'_640.tif');
A1=tiffread(str,[1 VH.ImageNumber]);
set(VH.text10,'string','Generating ROI ...')
mywaitbar(0,VH.axes2,'');
ROIfit1=single(zeros(2*fitl+1,2*fitl+1,m));
for i=1:m
    n=length(VH.DV(i,1).trackInfo(:,1));
    cx=round(VH.DV(i,1).trackInfo(1,3));
    cy=round(VH.DV(i,1).trackInfo(1,2));
    up=cx-fitl; 
    bottom=cx+fitl; 
    left=cy-fitl;
    right=cy+fitl;   
    if up>=1 && left>=1 && bottom<=row  && right<=column  && n<VH.ImageNumber             
        frame=VH.DV(i,1).trackInfo(:,1); 
        I=double(A1(up:bottom,left:right,frame));
        ROI=sum(I,3);
%         ROI=ROI-min(ROI(:));  
        ROIfit1(:,:,i)=single(ROI);            
    end 
    if (floor(i/m*100))>(floor((i-1)/m*100))
    mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
end

set(VH.text10,'string','GPU fitting ...')
N=floor(m/10000);
for i=1:N+1
    st=(i-1)*10000+1;
    et=i*10000;
    et=min(et,m);
    ROI=ROIfit1(:,:,st:et);%*factor/gain;
    tic
    [X Y N BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(single(ROI),sigma,20,2);
    toc
%     num = length(X);
%     M3 = zeros(num,1);
%     for j=1:num
%         ccx = X(j,1)-fitl;
%         ccy = Y(j,1)-fitl;
%         sx = median(S);
%         ROI1 = (exp(-0.5*(Xd-ccx).^2./(sx^2)-0.5*(Yd-ccy).^2./(sx^2)));
%         Dist = ((Xd-ccx).^2+(Yd-ccy).^2).^(3/2);
%         ROI1 = ROI1.*Dist;
%         ROI1 = ROI1/sum(ROI1(:));
%         ROI1 = ROI1.*single(ROI(:,:,j));
%         ROI1 = ROI1.*(ROI1>0);
%         M3(j,1) = sum(ROI1(:));
%     end
%     fitInfo2(st:et,3)=M3;
    fitInfo1(st:et,4)=N;
    CX(st:et,4)=X+1;
    CY(st:et,4)=Y+1;
end
%**************************************************************************
ROIfit1=[];
ROI = [];
ROI1 = [];

M0 = fitInfo1;
% meanI = mean(M0);
% for i=1:4
%     M0(:,i) = M0(:,i)*max(meanI)/meanI(i);
% end

% M0(:,1)=M0(:,1)*1.61;
% M0(:,2)=M0(:,2);
% M0(:,3)=M0(:,3)*1.26;
% M0(:,4)=M0(:,4)*2.79;

fitInfo1 = [];
fitInfo1(:,1) = M0(:,1)./(M0(:,1)+M0(:,4));
fitInfo1(:,4) = M0(:,4)./(M0(:,1)+M0(:,4));
fitInfo1(:,2) = M0(:,2)./(M0(:,2)+M0(:,3));
fitInfo1(:,3) = M0(:,3)./(M0(:,2)+M0(:,3));
% fitInfo1(:,1) = (fitInfo1(:,1)-min(fitInfo1(:,1)))/(max(fitInfo1(:,1)-min(fitInfo1(:,1))));
% fitInfo1(:,2) = (fitInfo1(:,2)-min(fitInfo1(:,2)))/(max(fitInfo1(:,2)-min(fitInfo1(:,2))));
% fitInfo1(:,3) = (fitInfo1(:,3)-min(fitInfo1(:,3)))/(max(fitInfo1(:,3)-min(fitInfo1(:,3))));
% fitInfo1(:,4) = (fitInfo1(:,4)-min(fitInfo1(:,4)))/(max(fitInfo1(:,4)-min(fitInfo1(:,4))));
zpos1 = MyCalPahse(fitInfo1(:,1),fitInfo1(:,4),fitInfo1(:,2),fitInfo1(:,3),0,60/180*pi)*300;

% M3 = fitInfo2;
% fitInfo2 = [];
% fitInfo2(:,1) = M3(:,1)./(M3(:,1)+M3(:,4));
% fitInfo2(:,4) = M3(:,4)./(M3(:,1)+M3(:,4));
% fitInfo2(:,2) = M3(:,2)./(M3(:,2)+M3(:,3));
% fitInfo2(:,3) = M3(:,3)./(M3(:,2)+M3(:,3));
% zpos2 = MyCalPahse(fitInfo2(:,1),fitInfo2(:,4),fitInfo2(:,2),fitInfo2(:,3),0,70/180*pi)*300;

zpos=zpos1;
% zpos = zeros(m,1);
% IX = zpos2-zpos1;
% zpos = zpos1.*(abs(IX)<50);
% for i=1:m
%     if IX(i,1)==0
%         zpos(i,1)=0;
%     end
%     if  (abs(IX(i,1))<50 || abs(abs(IX(i,1))-300)<50)
%         zpos(i,1)=zpos1(i,1);
%     end
% %     if IX(i,1)>0 && (abs(IX(i,1))<50 || abs(abs(IX(i,1))-300)<50)
% %         zpos(i,1)=zpos1(i,1)-300;
% %     end
% end
% M0_fit(:,1) = (M0(:,1)-M0(:,4))./(M0(:,1)+M0(:,4));
% M0_fit(:,2) = (M0(:,2)-M0(:,3))./(M0(:,2)+M0(:,3));
% M0_fit = (M0_fit+1)/2;
% zpos = sim(P.net,M0_fit');
M0 = [];
M3 = [];
fitInfo1 = [];
fitInfo2 = [];
%**************************************************************************

% for  i=1:L
%     
% end
CX(:,5)=fitInfo(:,1);
CY(:,5)=fitInfo(:,2);
fitInfo(:,1)=median(CX,2);
fitInfo(:,2)=median(CY,2);
VH.fitInfo_old=fitInfo;

% VH.ROI_old=ROIfit;
mm=1;
VH.ROI=uint16(zeros(fitsize,fitsize,m));
% VH.ROIfit=cell(m,1);
VH.trackInfo=cell(m,1);
VH.fitInfo=zeros(m,9);
% [X Y]= meshgrid(1:2*fitl2+1,1:2*fitl2+1);
set(VH.text10,'string','Generating information ...')
mywaitbar(0,VH.axes2,'');
for i=1:m
    cx=VH.DV(i,1).trackInfo(1,3);
    cy=VH.DV(i,1).trackInfo(1,2);
    n=length(VH.DV(i,1).trackInfo(:,1));
    up=round(cx)-fitl; 
    bottom=round(cx)+fitl; 
    left=round(cy)-fitl;
    right=round(cy)+fitl;       
    if  up>=1 && left>=1 && bottom<=row  && right<=column  && n<VH.ImageNumber && abs(fitInfo(i,1)-fitl-1)<1 && abs(fitInfo(i,2)-fitl-1)<1 ...
        && fitInfo(i,1)~=round(fitInfo(i,1)) && fitInfo(i,2)~=round(fitInfo(i,2)) && zpos(i,1)>0
        I=ROIfit(:,:,i);
        VH.ROI(:,:,mm)=I;
        if photonflag==1 || fitInfo(i,3)==0
            Nm=sum(I(:))*factor/gain;
        else
            Nm=fitInfo(i,3);
        end
        peak=Nm;
        
        if  fitInfo(i,5)>minSigma && fitInfo(i,5)<maxSigma
            sx=fitInfo(i,5);
        else
            I=I-(min(min(I))+max(max(I)))/2;
            I=I>0;
            se=strel('square',1);
            I=imclose(I,se);
            N=nnz(I);
            sx=sqrt(N)/2;
        end            

        if gap>0
            last_frame=VH.DV(i,1).trackInfo(n,1);
            first_frame=VH.DV(i,1).trackInfo(1,1);
            rcx=round(cx);
            rcy=round(cy);
            if first_frame==1
                I=double(VH.A(rcx-fitl:rcx+fitl,rcy-fitl:rcy+fitl,last_frame+1));
                sd=std2(I);
            elseif last_frame==VH.ImageNumber
                I=double(VH.A(rcx-fitl:rcx+fitl,rcy-fitl:rcy+fitl,first_frame-1)); 
                sd=std2(I);
            else
                I1=double(VH.A(rcx-fitl:rcx+fitl,rcy-fitl:rcy+fitl,first_frame-1)); 
                I2=double(VH.A(rcx-fitl:rcx+fitl,rcy-fitl:rcy+fitl,last_frame+1));
                sd=(std2(I1)+std2(I2))/2;
            end
        else
            sd=200; 
        end        
        b=sd*factor/gain;  

         if abs(fitInfo(i,1)-fitl-1-(cx-round(cx)))<0.5 && abs(fitInfo(i,2)-fitl-1-(cy-round(cy)))<0.5 && abs(fitInfo(i,1)-fitl-1)>0 && abs(fitInfo(i,2)-fitl-1)>0 && abs(fitInfo(i,1)-fitl-1)<1 && abs(fitInfo(i,2)-fitl-1)<1
            cx=round(cx)+fitInfo(i,1)-fitl-1;
            cy=round(cy)+fitInfo(i,2)-fitl-1;    
%             ccx=fitInfo(i,1)+fitl2-fitl;
%             ccy=fitInfo(i,2)+fitl2-fitl;
%         else
%             ccx=fitl2+1+cx-round(cx);
%             ccy=fitl2+1+cy-round(cy);
        end
%         VH.ROIfit{mm}=Nm*(exp(-0.5*(X-ccx).^2./(sx^2)-0.5*(Y-ccy).^2./(sx^2)));
        s=sx*a;                  
        if sx==0
            sx=1;
        end
        deltaR=sqrt((s^2+a^2/12)/Nm+8*pi*s^4*b^2/a^2/Nm^2);
        zp = zpos(i,1);
        VH.fitInfo(mm,:)=[cx cy peak sx sx deltaR Nm b zp];
        VH.trackInfo{mm}=VH.DV(i,1).trackInfo;
        mm=mm+1;
    end
    if (floor(i/m*100))>(floor((i-1)/m*100))
    mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
end