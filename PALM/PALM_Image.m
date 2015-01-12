clear
clc
fitl2=3;
sx=1;
[X Y]= meshgrid(1:2*fitl2+1,1:2*fitl2+1);
ROI=1000*(exp(-0.5*(X-fitl2-1).^2./(sx^2)-0.5*(Y-fitl2-1).^2./(sx^2)));
ROI=ROI/sum(ROI(:))*1000*300/10.5;
A=tiffread('F:\rapidStorm\Test.tif');
% A=zeros(128,128);
% for i=1:128
%     for j=1:128
%         dist=sqrt((i-64)^2+(j-64)^2);
%         if (dist<=10 && dist>=8) || (dist<=14 && dist>=12) || (dist<=18 && dist>=16) || (dist<=22 && dist>=20) || (dist<=26 && dist>=24) || (dist<=30 && dist>=28) || (dist<=34 && dist>=32) ...
%                 || (dist<=38 && dist>=36) || (dist<=42 && dist>=40)  || (dist<=46 && dist>=44) || (dist<=50 && dist>=48) 
%             A(i,j)=1;
%         end
%     end
% end
[row,col]=size(A);
[y,x]=find(A>0);
X=[x,y];
n=length(X);
for i=1:1
    Xnew((i-1)*n+1:i*n,:)=X;
end

X=Xnew;
n=length(X);
I=uint16(zeros(row,col,200));
h=waitbar(0,'Please wait ...');
fn=round(n/200);
V=[];
for i=1:200
    J=double(zeros(row,col));
    for j=1:fn
        N=round(rand(1)*n);
        N=max(N,1);
        cx=X(N,1);
        cy=X(N,2);
        flag=0;  
        if i>1
            CX=V{i-1}(:,1);
            CY=V{i-1}(:,2);
            dist=(cx-CX).^2+(cy-CY).^2;
            if any(dist<49)
                flag=1;
            end
        end
        if j>1
            CX1=V{i}(1:j-1,1);
            CY1=V{i}(1:j-1,2);
            dist=(cx-CX1).^2+(cy-CY1).^2;
            if any(dist<49)
                flag=1;
            end
        end
        while flag==1
            N=round(rand(1)*n);
            N=max(N,1);
            cx=X(N,1);
            cy=X(N,2);
            flag=0;
            if i>1
                CX=V{i-1}(:,1);
                CY=V{i-1}(:,2);
                dist=(cx-CX).^2+(cy-CY).^2;
                if any(dist<49)
                    flag=1;
                end
            end
            if j>1
                CX1=V{i}(1:j-1,1);
                CY1=V{i}(1:j-1,2);
                dist=(cx-CX1).^2+(cy-CY1).^2;
                if any(dist<49)
                    flag=1;
                end
            end
        end        
        J(cy-fitl2:cy+fitl2,cx-fitl2:cx+fitl2)=J(cy-fitl2:cy+fitl2,cx-fitl2:cx+fitl2)+ROI;  
        V{i}(j,1)=cx;
        V{i}(j,2)=cy;
    end
    J=J+0.0;
%     J=imnoise(J,'gaussian',0.000,0.0001);
%     J=imnoise(J,'poisson');
    I(:,:,i)=J;
    waitbar(i/100)
end
close(h)
dipstart
clc
I=noise(I,'poisson',0.035);
I=dip_array(I);
I=I+150;
I=noise(I,'gaussian',20);
I=dip_array(I);
tiffwrite(I,'F:\rapidStorm\test1.tif');

% J=double(zeros(row,col,200));
% h=waitbar(0,'Please wait ...');
% for i=1:length(X)
%     waitbar(i/length(X))
%     sf=round(200*rand(1));
%     sf=max(sf,1);
%     l=round(10*(0.2+0.5*randn(1)));
%     l=max(l,1);
%     ef=sf+l-1;
%     ef=min(ef,200);
%     cx=X(i,1);
%     cy=X(i,2);
%     for j=sf:ef
%         J(cy-3:cy+3,cx-3:cx+3,j)=J(cy-3:cy+3,cx-3:cx+3,j)+ROI;      
%     end
% end
% close(h)
% h=waitbar(0,'Please wait ...');
% for i=1:200
%     waitbar(i/200)
%     K=J(:,:,i)+0.1;
%     K=imnoise(K,'gaussian',0.01,0.01);
%     K=imnoise(K,'poisson');
%     I(:,:,i)=uint16(K*1000);
% end
% close(h)
% tiffwrite(I);
