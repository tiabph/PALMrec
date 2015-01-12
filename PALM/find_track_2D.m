clc
n=length(handles.DV);
[row,col]=size(handles.A(:,:,1));
m=1;
R=[];
B=[];
st=1;
VG=[];
VM=[];
VZ=[];
PN=[];
fitl=3;
[Xd Yd]= meshgrid(-3:3,-3:3);
xdata=[Xd Yd];
% pathname='E:\test\alxea647\';
            
for i=1:n
    l=length(handles.DV(i,1).trackInfo(:,1));
    cy=round(handles.DV(i,1).trackInfo(1,3));
    cx=round(handles.DV(i,1).trackInfo(1,2));
    up=round(cy)-fitl; 
    bottom=round(cy)+fitl; 
    left=round(cx)-fitl;
    right=round(cx)+fitl;       
    ROI=[];
    if l>5 && l<30 && up>=1 && left>=1 && bottom<=row  && right<=col && all(handles.DV(i,1).trackInfo(:,4)>200) %&& mean(handles.DV(i,1).trackInfo(:,4))<1000
        for j=1:l
            f=handles.DV(i,1).trackInfo(j,1);
            ROI(:,:,j)=double(handles.A(cy-3:cy+3,cx-3:cx+3,f));
        end
        
        C=[];
        C=handles.DV(i,1).trackInfo(:,2:3);
        C(:,1)=C(:,1)-mean(C(:,1));
        C(:,2)=C(:,2)-mean(C(:,2));
        
        X=[];
        Y=[];
        [X Y N1 BG S CRLBx CRLBy CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(single(ROI),1.2,20,1);
        B(m,1)=mean(X)-3;
        B(m,2)=mean(Y)-3;
        X=X-mean(X);
        Y=Y-mean(Y);
        
        if  all(abs(C(:))<0.5) && all(abs(X)<0.5) && all(abs(Y)<0.5)
            et=st+l-1;
            VM(st:et,1)=X;
            VM(st:et,2)=Y;
            st=st+l;
            m=m+1
        end            
    end
end

VM=VM*107;
figure;
plot(VM(:,1),VM(:,2),'b.');
std(VM)