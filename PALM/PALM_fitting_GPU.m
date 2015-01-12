function VH=PALM_fitting_GPU(VH)

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
meanBg=VH.parameter.reconstruction.meanBg;

row=VH.row;
column=VH.column;
m=length(VH.DV);
% fitl=round(3*sigma);

set(VH.text10,'string','Generating ROI ...')
mywaitbar(0,VH.axes2,'');
ROIfit=single(zeros(2*fitl+1,2*fitl+1,m));
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
    if up>=1 && left>=1 && bottom<=row  && right<=column  && n<VH.ImageNumber && abs(fitInfo(i,1)-fitl-1)<1 && abs(fitInfo(i,2)-fitl-1)<1 && fitInfo(i,1)~=round(fitInfo(i,1)) && fitInfo(i,2)~=round(fitInfo(i,2))
        I=ROIfit(:,:,i);
        VH.ROI(:,:,mm)=I;
        if photonflag==0           
            Nm=sum(I(:))*factor/gain;
        else
            Nm=fitInfo(i,3);
        end
        peak=fitInfo(i,3);
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
        
        last_frame=VH.DV(i,1).trackInfo(n,1);
        first_frame=VH.DV(i,1).trackInfo(1,1);
        if gap>0
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
            if meanBg>0
                b=meanBg;
            else
                b=5;
            end
            sd=b*gain/factor;
        end
        b=sd*factor/gain;  

        if abs(fitInfo(i,1)-fitl-1)<1 && abs(fitInfo(i,2)-fitl-1)<1 && abs(fitInfo(i,1)-fitl-1)>0 && abs(fitInfo(i,2)-fitl-1)>0 ...
            && abs(fitInfo(i,1)-fitl-1-(cx-round(cx)))<0.5 && abs(fitInfo(i,2)-fitl-1-(cy-round(cy)))<0.5
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
        deltaR=sqrt((s^2+a^2/12)/Nm+8*pi*s^4*b^2/a^2/Nm^2);
        VH.fitInfo(mm,:)=[cx cy peak sx sx deltaR Nm b first_frame];
        VH.trackInfo{mm}=VH.DV(i,1).trackInfo;
        mm=mm+1;
    end
    if (floor(i/m*100))>(floor((i-1)/m*100))
    mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
end