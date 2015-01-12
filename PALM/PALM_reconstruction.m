function PALM_reconstruction(VH)

amp1=VH.parameter.reconstruction.amp;
maxError=VH.parameter.reconstruction.maxError;
minPhoton=VH.parameter.reconstruction.minPhoton;
maxPhoton=VH.parameter.reconstruction.maxPhoton;
maxBg=VH.parameter.reconstruction.maxBg;
meanBg=VH.parameter.reconstruction.meanBg;
minSigma=VH.parameter.reconstruction.minSigma;
maxSigma=VH.parameter.reconstruction.maxSigma;
meanSigma=VH.parameter.reconstruction.meanSigma;
maxLength=VH.parameter.reconstruction.maxLength;
drift=VH.parameter.driftCorrection;

if VH.parameter.prompt_flag
    VH.parameter.batch=0;
end

if VH.parameter.PALM.fitting_flag==0;
    set(VH.text10,'string','Loading data ...')   
    pause(0.1)
    if VH.parameter.fitting.flag==3  
        fstr=strcat(VH.pathname,VH.filebase,'_Fitting_ANN.mat');
    else
        fstr=strcat(VH.pathname,VH.filebase,'_Fitting.mat');
    end
    if exist(fstr,'file')
        P=importdata(fstr);
        VH.fitInfo=P.fitInfo;
        if VH.parameter.fitting.flag==3  
            VH.fitInfo_old=P.fitInfo_old;
        end
        VH.trackInfo=P.trackInfo;
        VH.ROI=P.ROI;
        fitl=(P.parameter.fitting.fitsize-1)/2;
        fitl2=fitl+2;
        a=P.parameter.fitting.pixelsize;
        number=P.ImageNumber;
        if ~isfield(VH,'row')
            VH.row=P.row;
            VH.column=P.column;
        end
        P=[];
    end
else
    fitl=(VH.parameter.fitting.fitsize-1)/2;
    fitl2=fitl+2;
    a=VH.parameter.fitting.pixelsize;
    number=VH.ImageNumber;
end

if drift==1
    str=strcat(VH.pathname,VH.filebase,'_beads.mat');
    if exist(str,'file')
        P=importdata(str);
        B=P;
    else
        str=strcat(VH.pathname,VH.filebase,'_beads.tif');
        if exist(str,'file')
            I=tiffread(str,[1 number]);
            set(VH.text10,'string','Fitting beads ...')
            pause(eps)
            n=length(I);
            wid=size(I,1);
            [Xd,Yd] = meshgrid(1:wid,1:wid);
            xdata=[Xd,Yd];
            B=zeros(n,2);
            mywaitbar(0,VH.axes2,'');
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
                [yy xx]=weight_centrid(Fd,wid);
                if abs(yy-cy)>1 || abs(xx-cx)>1
                    yy=cy;
                    xx=cx;
                end
                if abs(B(i,1)-yy)>2 || abs(B(i,2)-xx)>2
                    B(i,1)=yy;
                    B(i,2)=xx;
                end
                if (floor(i/n*100))>(floor((i-1)/n*100))
                    mywaitbar(i/n,VH.axes2,[num2str(floor(i/n*100)),'%']);
                end
            end

            for i=2:length(B)-1
                if abs(B(i,1)-B(i-1,1))>3 && abs(B(i-1,1)-B(i+1,1))<0.5
                    B(i,1)=(B(i-1,1)+B(i+1,1))/2;
                end
                if abs(B(i,2)-B(i-1,2))>3 && abs(B(i-1,2)-B(i+1,2))<0.5
                    B(i,2)=(B(i-1,2)+B(i+1,2))/2;
                end
            end

            str=strcat(VH.pathname,'beads.tif');
            if exist(str,'file')
                I=tiffread(str);
                Fd=double(I);
                backg=min(Fd(:));
                peak=max(Fd(:));
                [cy cx]=find(Fd==peak);
                cy=mean(cy);
                cx=mean(cx);
                xp=[cx,cy,2,2,backg,peak];
                options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',1000,'TolFun',0.01,'Algorithm','levenberg-marquardt');
                [lp,resnorm,residual,exitflag]=lsqcurvefit(@Gaussian_PALM1,xp,xdata,Fd,[],[],options);
                x=lp(1);
                y=lp(2);
                [yy xx]=weight_centrid(Fd,wid);
                if abs(yy-cy)>1 || abs(xx-cx)>1
                    yy=cy;
                    xx=cx;
                end
                if abs(y-yy)>2 || abs(x-xx)>2
                    y=yy;
                    x=xx;
                end
                B(:,1)=B(:,1)-y;
                B(:,2)=B(:,2)-x;
            else
%                 B(:,1)=B(:,1)-B(end,1);
%                 B(:,2)=B(:,2)-B(end,2);
                B(:,1)=B(:,1)-B(1,1);
                B(:,2)=B(:,2)-B(1,2);
            end
%             B(:,1)=B(:,1)+0.2723;
%             B(:,2)=B(:,2)+0.5450;            
            B(:,1)=smooth(B(:,1),50,'rlowess');
            B(:,2)=smooth(B(:,2),50,'rlowess');
            if ~VH.parameter.batch
                figure;plot(B(:,1),'b');hold on;plot(B(:,2),'r');
            end
            str=strcat(VH.pathname,VH.filebase,'_beads.mat');
            save(str,'B');
        else
            drift=0;
        end       
    end
end

set(VH.text10,'string','Recontructing image ...')
row=VH.row;
column=VH.column;
amp2=8;
row_new1=row*amp1;
column_new1=column*amp1;
m=sum(VH.fitInfo(:,4)>0);
% if meanSigmaflag>0
% sx=zeros(m,1);
% mywaitbar(0,VH.axes2,'');
% for i=1:m   
%     sx(i,1)=VH.fitInfo(i,4);
%     if (floor(i/m*100))>(floor((i-1)/m*100))
%         mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
%     end
% end
% meanSigma=mean(sx);
% else
%     meanSigma=VH.parameter.fitting.sigma;
% end
[X Y]= meshgrid(1:2*fitl2+1,1:2*fitl2+1);
% ROI= 1000*(exp(-0.5*(X-fitl2-1).^2./(sx^2)-0.5*(Y-fitl2-1).^2./(sx^2)));
row_new2=row*amp2;
column_new2=column*amp2;
Reconstruct_I=single(zeros(row_new1,column_new1));
Reconstruct_I1=single(zeros(row_new2,column_new2));
Reconstruct_I2=single(zeros(row_new1,column_new1));
Reconstruct_I3=single(zeros(row_new2,column_new2));
Reconstruct_I6=single(zeros(row_new1,column_new1));
Reconstruct_I_3D=single(zeros(row_new2,column_new2,150));
% Reconstruct_I_3D2=single(zeros(row_new2,column_new2,150));
% if VH.parameter.fitting.flag==3  
%     Reconstruct_I4=double(zeros(row_new1,column_new1,3));
%     Reconstruct_I7=double(zeros(row_new1,column_new1,3));
% %     Reconstruct_I5=double(zeros(row_new1,column_new1));
%     map=colormap(jet(101));
% end
peak=zeros(m,1);
sx=zeros(m,1);
deltaR=zeros(m,1);
photonN=zeros(m,1);
backN=zeros(m,1);
V=zeros(m,3);
n=0;

if maxLength==0
    maxLength=VH.ImageNumber;
end

% str='z_hist.mat';
% P=importdata(str);
mywaitbar(0,VH.axes2,'');
for i=1:m
    l=length(VH.trackInfo{i,1}(:,1));
    t=VH.trackInfo{i,1}(1,1);
     if  l<=maxLength && VH.fitInfo(i,7)>minPhoton && VH.fitInfo(i,7)<maxPhoton && VH.fitInfo(i,6)<maxError && VH.fitInfo(i,8)<maxBg && VH.fitInfo(i,4)>minSigma && VH.fitInfo(i,4)<maxSigma
         t=VH.trackInfo{i,1}(1,1);
         n=n+1;
         peak(n,1)=VH.fitInfo(i,3);
         sx(n,1)=VH.fitInfo(i,4);
         s=sx(n,1)*a;
         if meanSigma>0
             s=meanSigma*a;
         end         
         backN(n,1)=VH.fitInfo(i,8);
         if meanBg>0
             b=meanBg;
         else
             b=VH.fitInfo(i,8);
         end
         Nm=VH.fitInfo(i,7);
         deltaR(n,1)=sqrt((s^2+a^2/12)/Nm+8*pi*s^4*b^2/a^2/Nm^2);
         photonN(n,1)=VH.fitInfo(i,7);
         if drift==1            
            VH.fitInfo(i,1)=VH.fitInfo(i,1)-B(t,1);
            VH.fitInfo(i,2)=VH.fitInfo(i,2)-B(t,2);
         end
         ccx=round(VH.fitInfo(i,1)*amp1);
         ccy=round(VH.fitInfo(i,2)*amp1);    
         if ccx-fitl2>0 && ccy-fitl2>0 && ccx+fitl2<=row_new1 && ccy+fitl2<=column_new1
             ccx1=fitl2+1+VH.fitInfo(i,1)*amp1-ccx;
             ccy1=fitl2+1+VH.fitInfo(i,2)*amp1-ccy;
             ROI=1000*(exp(-0.5*(X-ccx1).^2./(1^2)-0.5*(Y-ccy1).^2./(1^2)));
             ROI1=Nm*(exp(-0.5*(X-ccx1).^2./(sx(n,1)^2)-0.5*(Y-ccy1).^2./(sx(n,1)^2)));
             Reconstruct_I(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2)=Reconstruct_I(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2)+ROI1;
             Reconstruct_I2(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2)=Reconstruct_I2(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2)+ROI;
             Reconstruct_I6(ccx,ccy)=Reconstruct_I6(ccx,ccy)+1;
%             if VH.parameter.fitting.flag==3  
%                 if VH.fitInfo_old(i,7)>0 % VH.fitInfo_old(i,8)>=0 && VH.fitInfo_old(i,8)<=1 %
%                     delta=round(VH.fitInfo_old(i,7));
%                     mobility=ceil(sind(delta/2)*100);
%                     rgb=map(mobility+1,:);
% %                     asymetry=ceil(VH.fitInfo_old(i,8)*100);
% %                     rgb=map(delta+1,:);
% %                     rgb=map(asymetry,:);
%                     ROI=1*(exp(-0.5*(X-ccx1).^2./(meanSigma^2)-0.5*(Y-ccy1).^2./(meanSigma^2)));
%                     Reconstruct_I4(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2,1)=Reconstruct_I4(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2,1)+ROI*(rgb(1));
%                     Reconstruct_I4(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2,2)=Reconstruct_I4(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2,2)+ROI*(rgb(2));
%                     Reconstruct_I4(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2,3)=Reconstruct_I4(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2,3)+ROI*(rgb(3));
%                     Reconstruct_I7(ccx,ccy,1)=Reconstruct_I7(ccx,ccy,1)+rgb(1);
%                     Reconstruct_I7(ccx,ccy,2)=Reconstruct_I7(ccx,ccy,2)+rgb(2);
%                     Reconstruct_I7(ccx,ccy,3)=Reconstruct_I7(ccx,ccy,3)+rgb(3);
% %                 else
% %                     Reconstruct_I5(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2)=Reconstruct_I5(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2)+ROI;
%                 end
%             end
         end
         ccx=round(VH.fitInfo(i,1)*amp2);
         ccy=round(VH.fitInfo(i,2)*amp2);
         V(n,1)=VH.fitInfo(i,1);
         V(n,2)=VH.fitInfo(i,2);
         V(n,3)=t;
         V(n,4)=Nm;
         if ccx-fitl>0 && ccy-fitl>0  && ccx+fitl<=row_new2 && ccy+fitl<=column_new2
             Reconstruct_I1(ccx-fitl:ccx+fitl,ccy-fitl:ccy+fitl)=Reconstruct_I1(ccx-fitl:ccx+fitl,ccy-fitl:ccy+fitl)+double(VH.ROI(:,:,i));
             Reconstruct_I3(ccx,ccy)=Reconstruct_I3(ccx,ccy)+1;%double(VH.ROI(fitl2+1,fitl2+1,i));
%              Reconstruct_I6(ccx,ccy)=Reconstruct_I6(ccx,ccy)+1;
             if VH.parameter.fitting.flag==4      
                 if drift==1
                    VH.fitInfo(i,9)=VH.fitInfo(i,9)-B(t,3); 
                 end
                 z=VH.fitInfo(i,9);
                 if z>300
                     z=z-300;
                 end
                 if z<0
                     z=z+300;
                 end
                 V(n,5)=z;
%                  tt=ceil(t/100);
%                  VH.fitInfo(i,9)=VH.fitInfo(i,9)-P(tt,1);
                 if z>0 && z<=300      
                     z=round(z/2);
                     if z<1
                         z=1;
                     end    
                         Reconstruct_I_3D(ccx,ccy,z)=Reconstruct_I_3D(ccx,ccy,z)+100;
%                          ccx1=fitl2+1+VH.fitInfo(i,1)*amp2-ccx;
%                          ccy1=fitl2+1+VH.fitInfo(i,2)*amp2-ccy;
%                          ROI=1000*(exp(-0.5*(X-ccx1).^2./(meanSigma^2)-0.5*(Y-ccy1).^2./(meanSigma^2)));     
%                          Reconstruct_I_3D2(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2,z)=Reconstruct_I_3D2(ccx-fitl2:ccx+fitl2,ccy-fitl2:ccy+fitl2,z)+ROI;
%                      end
                 end
             end
         end
     end
     if (floor(i/m*100))>(floor((i-1)/m*100))
         mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
     end
end

% peak=peak(peak>0);
V=V(V(:,1)>0,:);
photonN=photonN(photonN>0);
deltaR=deltaR(deltaR>0);
sx=sx(sx>0);
backN=backN(backN>0);
A=[];
A(1)=mean(photonN);
A(2)=median(photonN);
A(3)=mean(deltaR);
A(4)=median(deltaR);
A(5)=mean(sx);
A(6)=mean(backN);
disp(['mean photons:',num2str(A(1))])
disp(['median photons:',num2str(A(2))])
disp(['mean position error:',num2str(A(3))])
disp(['median position error:',num2str(A(4))])
disp(['mean sigma:',num2str(A(5))]);
disp(['mean background photons:',num2str(A(6))])
file=strcat(VH.pathname,VH.filebase,'.txt');
fid=fopen(file,'w');
fprintf(fid,'%f\r\n',A);
fclose(fid);

file=strcat(VH.pathname,VH.filebase,'_position.txt');
dlmwrite(file, V, 'delimiter','\t','precision', 6);

if ~VH.parameter.batch
    figure;
    subplot(2,2,1);hist(photonN,100);title('total photons');
    subplot(2,2,2);hist(deltaR,100);title('position error');
    subplot(2,2,3);hist(sx,20);title('standard deviation');
    subplot(2,2,4);hist(backN,100);title('background photons');
end

warning off;
maxp=max(max(Reconstruct_I));
if maxp>=65535
minp=min(min(Reconstruct_I));
Reconstruct_I=(Reconstruct_I-minp)/(maxp-minp)*65535;
end
Reconstruct_I=uint16(Reconstruct_I);

maxp=max(max(Reconstruct_I1));
if maxp>=65535
minp=min(min(Reconstruct_I1));
Reconstruct_I1=(Reconstruct_I1-minp)/(maxp-minp)*65535;
end
Reconstruct_I1=uint16(Reconstruct_I1);

maxp=max(max(Reconstruct_I2));
if maxp>=65535
minp=min(min(Reconstruct_I2));
Reconstruct_I2=(Reconstruct_I2-minp)/(maxp-minp)*65535;
end
Reconstruct_I2=uint16(Reconstruct_I2);

maxp=max(max(Reconstruct_I3));
if maxp>=65535
minp=min(min(Reconstruct_I3));
Reconstruct_I3=(Reconstruct_I3-minp)/(maxp-minp)*65535;
end
Reconstruct_I3=uint16(Reconstruct_I3);

filename=strcat(VH.pathname,VH.filebase);
% if drift==0
%     FileStr=strcat(filename,'_reconstruct_fit','_amp',int2str(amp1),'.tif');
% else
%     FileStr=strcat(filename,'_reconstruct_fit','_amp',int2str(amp1),'_DC.tif');
% end
% imwrite(Reconstruct_I,FileStr,'tif','Compression','none','WriteMode','overwrite')
if drift==0
    FileStr=strcat(filename,'_reconstruct_locate','_amp',int2str(amp1),'.tif');
else
    FileStr=strcat(filename,'_reconstruct_locate','_amp',int2str(amp1),'_DC.tif');
end
imwrite(Reconstruct_I2,FileStr,'tif','Compression','none','WriteMode','overwrite')
FileStr=strcat(filename,'_reconstruct_',int2str(amp2),'.tif');
imwrite(Reconstruct_I1,FileStr,'tif','Compression','none','WriteMode','overwrite')
FileStr=strcat(filename,'_reconstruct_dot_',int2str(amp2),'.tif');
imwrite(Reconstruct_I3,FileStr,'tif','Compression','none','WriteMode','overwrite')
FileStr=strcat(filename,'_reconstruct_dot_',int2str(amp1),'.tif');
imwrite(uint16(Reconstruct_I6),FileStr,'tif','Compression','none','WriteMode','overwrite')

if VH.parameter.fitting.flag==4  
    FileStr=strcat(filename,'_reconstruct_dot','_amp',int2str(amp2),'_3D.tif');
    tiffwrite(Reconstruct_I_3D,FileStr);
%     FileStr=strcat(filename,'_reconstruct_locate','_amp',int2str(amp2),'_3D.tif');
%     tiffwrite(Reconstruct_I_3D2,FileStr);
end

% Reconstruct_I4_old=Reconstruct_I4;
% if VH.parameter.fitting.flag==3  
%     maxp=max(max(Reconstruct_I4(:)));
%     minp=min(min(Reconstruct_I4(:)));
%     Reconstruct_I4=(Reconstruct_I4-minp)/(maxp-minp);
%     
%     maxp=max(max(Reconstruct_I7(:)));
%     minp=min(min(Reconstruct_I7(:)));
%     Reconstruct_I7=(Reconstruct_I7-minp)/(maxp-minp);
%     
% %     maxp=max(max(Reconstruct_I5));
% %     if maxp>65535
% %         minp=min(min(Reconstruct_I5));
% %         Reconstruct_I5=(Reconstruct_I5-minp)/(maxp-minp)*65535;
% %     end
% %     Reconstruct_I5=uint16(Reconstruct_I5);
%     if drift==0
%         FileStr1=strcat(filename,'_reconstruct_ANN_restricted','_amp',int2str(amp1),'.tif');
% %         FileStr2=strcat(filename,'_reconstruct_ANN_free','_amp',int2str(amp1),'.tif');
%     else
%         FileStr1=strcat(filename,'_reconstruct_ANN_restricted','_amp',int2str(amp1),'_DC.tif');
% %         FileStr2=strcat(filename,'_reconstruct_ANN_free','_amp',int2str(amp1),'_DC.tif');
%     end   
%     imwrite(Reconstruct_I4,FileStr1,'tif','Compression','none','WriteMode','overwrite')
% %     imwrite(Reconstruct_I5,FileStr2,'tif','Compression','none','WriteMode','overwrite')
%     FileStr=strcat(filename,'_reconstruct_dot_new_',int2str(amp1),'.tif');
%     imwrite(Reconstruct_I7,FileStr,'tif','Compression','none','WriteMode','overwrite')
% end

fstr=strcat(VH.pathname,VH.filebase,'_Reconstruction.mat');
set(VH.text10,'string','Saving ...')
VH=[];
pause(eps)
save(fstr);

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