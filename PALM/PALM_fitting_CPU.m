function VH=PALM_fitting_CPU(VH)

gain=VH.parameter.fitting.gain;
factor=VH.parameter.fitting.factor;
sigma=VH.parameter.fitting.sigma;
a=VH.parameter.fitting.pixelsize;
fitsize=VH.parameter.fitting.fitsize;
fitl=(VH.parameter.fitting.fitsize-1)/2;
photonflag=VH.parameter.fitting.photonflag;
gap=VH.parameter.linking.gap;
meanBg=VH.parameter.reconstruction.meanBg;

row=VH.row;
column=VH.column;
m=length(VH.DV);
mm=1;
set(VH.text10,'string','Gaussian fitting ...')
mywaitbar(0,VH.axes2,'');

VH.ROI=uint16(zeros(fitsize,fitsize,m));
% VH.ROIfit=[];
VH.fitInfo=zeros(m,9);
VH.fitInfo_old=zeros(m,6);
VH.trackInfo=cell(m,1);
% VH.ROI_old=[];
% VH.DVnew=[];
[X Y]= meshgrid(1:2*fitl+1,1:2*fitl+1);
% [X1 Y1]= meshgrid(1:2*fitl2+1,1:2*fitl2+1);
for i=1:m
%     ROI=zeros(fitl*2+1);
    n=length(VH.DV(i,1).trackInfo(:,1));
    cx=VH.DV(i,1).trackInfo(1,3);
    cy=VH.DV(i,1).trackInfo(1,2);
    up=round(cx)-fitl; 
    bottom=round(cx)+fitl; 
    left=round(cy)-fitl;
    right=round(cy)+fitl;       
        
    if up>=1 && left>=1 && bottom<=row  && right<=column  && n<VH.ImageNumber  
        last_frame=VH.DV(i,1).trackInfo(n,1);
        first_frame=VH.DV(i,1).trackInfo(1,1);
        rcx=round(cx);
        rcy=round(cy);
        rcx=max(fitl+1,rcx);
        rcy=max(fitl+1,rcy);
        rcx=min(rcx,row-fitl);
        rcy=min(rcy,column-fitl);
        if gap>0
            if last_frame<VH.ImageNumber
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
        bm=sd*factor/gain;    
        
        options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',100,'TolFun',0.01);
%         ROI1=[];
%         deltaR1=[];
        frame=VH.DV(i,1).trackInfo(:,1); 
        I=double(VH.A(up:bottom,left:right,frame));
        ROI=sum(I,3);
        peak=max(max(ROI))-min(min(ROI));
%         pI=max(ROI(:));
        bg=min(min(ROI));      
        xdata=[X Y];
%         VH.ROI_old(:,:,i)=ROI;
        initpar = double([fitl+1,fitl+1,sigma,sigma,bg,peak]);
        [lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian_PALM,initpar,xdata,ROI,[],[],options);
        VH.fitInfo_old(i,:)=lp;
        if (exitflag>=0 || exitflag==-4) && lp(6)>0 && lp(5)>0 && abs(lp(1)-fitl-1)<1 && abs(lp(2)-fitl-1)<1 && lp(3)>0 && lp(4)>0
            peak=lp(6);
%             I=ROI;
%             I=I-min(I(:));
%             I=I.*(I>0);
            I=ROI-min(ROI(:));
            VH.ROI(:,:,mm)=I;
            Nm=sum(I(:))*factor/gain;     
            if abs(lp(1)-fitl-1)<1 && abs(lp(2)-fitl-1)<1 && abs(lp(1)-fitl-1-(cx-rcx))<0.5 && abs(lp(2)-fitl-1-(cy-rcy))<0.5
                cx=rcx+lp(1)-fitl-1;
                cy=rcy+lp(2)-fitl-1;
                ccx=lp(1);%+fitl2-fitl;
                ccy=lp(2);%+fitl2-fitl;
            else
                ccx=fitl+1+cx-rcx;
                ccy=fitl+1+cy-rcy;
            end
            sx=lp(3);
            sy=lp(4);            
            if photonflag==1
                Nm=2*pi*lp(3)^2*peak*factor/gain;   
%                 ROIfit=peak*(exp(-0.5*(X-ccx).^2./(sx^2)-0.5*(Y-ccy).^2./(sx^2)));
%                 Nm=sum(ROIfit(:))*factor/gain;   
            end
            s=sx*a;                  
            deltaR=sqrt((s^2+a^2/12)/Nm+8*pi*s^4*bm^2/a^2/Nm^2);
            VH.fitInfo(mm,:)=[cx cy peak sx sy deltaR Nm bm first_frame];
            VH.trackInfo{mm}=VH.DV(i,1).trackInfo;
            mm=mm+1;
        end
    end
    if floor(i/m*100)>floor((i-1)/m*100)
    mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
end