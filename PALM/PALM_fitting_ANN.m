function VH=PALM_fitting_ANN(VH)

gain=VH.parameter.fitting.gain;
factor=VH.parameter.fitting.factor;
sigma=VH.parameter.fitting.sigma;
a=VH.parameter.fitting.pixelsize;
fitsize=VH.parameter.fitting.fitsize;
fitl=(VH.parameter.fitting.fitsize-1)/2;
photonflag=VH.parameter.fitting.photonflag;

row=VH.row;
column=VH.column;
m=length(VH.DV);

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
VH.ROI_old=ROIfit;

VH.maxI=zeros(1,m);
ROIfit=reshape(ROIfit,fitsize^2,m);
for i=1:m
    VH.maxI(1,i)=max(ROIfit(:,i));
    ROIfit(:,i)=ROIfit(:,i)/max(ROIfit(:,i));    
end
set(VH.text10,'string','ANN fitting ...')
pause(0.1);

type=3;
if type==1
    filestr='C:\Users\Administrator\Documents\MATLAB\ANN\net_fixed_noise_random.mat';
end
if type==3
    filestr='C:\Users\Administrator\Documents\MATLAB\ANN\net_restricted_noise_random_new.mat';
end
P=importdata(filestr);
outputs=ones(8,m)*nan;
outputs(1:2,:)=sim(P.net1,double(ROIfit));
outputs(3:4,:)=sim(P.net2,double(ROIfit));
outputs(3,:)=outputs(3,:).*VH.maxI*fitsize^2;
outputs(4,:)=outputs(4,:).*VH.maxI;
outputs(3,:)=(outputs(3,:)-outputs(4,:)*fitsize^2)*factor/gain;
if type==1
    out=sim(P.net3,double(ROIfit));    
    n=length(out(1,:));
    k=sqrt(out(1,:).^2+out(2,:).^2);
    k=[k
        k];
    out(1:2,:)=out(1:2,:)./k;
    for i=1:n
        if out(1,i)>=0 && out(2,i)>=0       % 0-90
            outputs(5,i)=asind(out(1,i));
        elseif out(1,i)>=0 && out(2,i)<0    % 90-180
            outputs(5,i)=180-asind(out(1,i));
        elseif out(1,i)<0 && out(2,i)>=0    % 270-360
            outputs(5,i)=360+asind(out(1,i));
        elseif out(1,i)<0 && out(2,i)<0     % 180-270
            outputs(5,i)=180-asind(out(1,i));
        end
    end
    outputs(6,:)=sim(P.net4,double(ROIfit));   
    outputs(6,outputs(6,:)>1)=1;
    outputs(6,outputs(6,:)<0)=0;
    outputs(6,:)=outputs(6,:)*90;
    outputs(7,:)=0;
    outputs(8,:)=sind(2*outputs(6,:));
end
if type==3
    out=sim(P.net3,double(ROIfit));    
    n=length(out(1,:));
    out(1,:)=out(1,:)./out(3,:);
    out(2,:)=out(2,:)./out(3,:);
    for i=1:n
    %     out(1,i)=min(max(out(1,i),-1),1);
    %     out(2,i)=min(max(out(2,i),-1),1);
        if out(1,i)<-1 || out(1,i)>1
            out(1,i)=nan;
        end
        if out(2,i)<-1 || out(2,i)>1
            out(2,i)=nan;
        end
    end
    for i=1:n
        if out(1,i)>=0 && out(2,i)>=0       % 0-90
            outputs(5,i)=asind(out(1,i));
        elseif out(1,i)>=0 && out(2,i)<0    % 90-180
            outputs(5,i)=180-asind(out(1,i));
        elseif out(1,i)<0 && out(2,i)>=0    % 270-360
            outputs(5,i)=360+asind(out(1,i));
        elseif out(1,i)<0 && out(2,i)<0     % 180-270
            outputs(5,i)=180-asind(out(1,i));
        end
    end
    out=sim(P.net4,double(ROIfit)); 
    outputs(8,:)=out(1,:);
    n=length(out(1,:));
    k=sqrt(out(1,:).^2+out(2,:).^2);
    for i=1:n
    %     k(1,i)=min(max(k(1,i),0),1);
        if k(1,i)==0 || k(1,i)>1
            k(1,i)=nan;
        end
    end
    outputs(7,:)=acosd(sqrt(2*k+1/4)-1/2);
    out(1,:)=out(1,:)./k;
    out(2,:)=out(2,:)./k;
    for i=1:n
        if out(1,i)>0 && out(2,i)>0       % 0-90
            outputs(6,i)=asind(out(1,i))/2;
        elseif out(1,i)>0 && out(2,i)<0    % 90-180
            outputs(6,i)=(180-asind(out(1,i)))/2;
        end
    end
end

VH.maxI=[];
P=[];
fitInfo=outputs';
VH.fitInfo_old=zeros(m,8);

ROIfit=VH.ROI_old;
VH.ROI_old=[];
VH.ROI=uint16(zeros(fitsize,fitsize,m));
% VH.ROIfit=cell(m,1);
VH.trackInfo=cell(m,1);
VH.fitInfo=zeros(m,8);
% [X Y]= meshgrid(1:2*fitl2+1,1:2*fitl2+1);
mm=1;
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
    if up>=1 && left>=1 && bottom<=row  && right<=column  && n<VH.ImageNumber
        I=ROIfit(:,:,i);
        VH.ROI(:,:,mm)=uint16(I);
        if photonflag==1
            Nm=sum(I(:))*factor/gain;
        else
            Nm=fitInfo(i,3);
        end
        peak=fitInfo(i,3);
        sx=sigma;                 
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
        b=sd*factor/gain;  
        if abs(fitInfo(i,1)-fitl-1)<1 && abs(fitInfo(i,2)-fitl-1)<1
            cx=cx+fitInfo(i,1)-fitl-1;
            cy=cy+fitInfo(i,2)-fitl-1;    
%             ccx=fitInfo(i,1)+fitl2-fitl;
%             ccy=fitInfo(i,2)+fitl2-fitl;
        else
%             ccx=fitl2+1+cx-round(cx);
%             ccy=fitl2+1+cy-round(cy);
        end
%         VH.ROIfit{mm}=Nm*(exp(-0.5*(X-ccx).^2./(sx^2)-0.5*(Y-ccy).^2./(sx^2)));
        s=sx*a;            
        deltaR=sqrt((s^2+a^2/12)/Nm+8*pi*s^4*b^2/a^2/Nm^2);
        VH.fitInfo(mm,:)=[cx cy peak sx sx deltaR Nm b];
        VH.fitInfo_old(mm,:)=fitInfo(i,:);
        VH.trackInfo{mm}=VH.DV(i,1).trackInfo;
        mm=mm+1;
    end
    if (floor(i/m*100))>(floor((i-1)/m*100))
        mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
end