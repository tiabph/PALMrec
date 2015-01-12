clear
clc
pathname='G:\zms\20130916_meos3.2F173S-lifeact\3\';
filebase='cell3';
N=4;    % stack number
r=2;    % distance (pixel)
t=0;    % time
for i=1:1
    i
    FI=[];
    m=0;
    
    str=strcat(pathname,filebase,'_l_',int2str(i),'_Fitting.mat');
    P=importdata(str);
    fitInfo=P.fitInfo;
    Vl=fitInfo(:,1:3);
    Vl(:,4)=fitInfo(:,9);
    P=[];
    fitInfo=[];
    str=strcat(pathname,filebase,'_r_',int2str(i),'_Fitting.mat');
    P=importdata(str);
    fitInfo=P.fitInfo;
    Vr=fitInfo(:,1:3);
    Vr(:,4)=fitInfo(:,9);
    P=[];
    fitInfo=[];   
    str=strcat(pathname,filebase,'_l_',int2str(i),'.tif');
    Il=tiffread(str);
    str=strcat(pathname,filebase,'_r_',int2str(i),'.tif');
    Ir=tiffread(str);
    
    Vr(:,5)=0;
    nl=sum(Vl(:,1)>0);
    nr=sum(Vr(:,1)>0);
    for j=1:nr    
        j
        IX=abs(Vl(:,4)-Vr(j,4))<=t;
        Vnew=Vl(IX,:);
        Dist=sqrt((Vnew(:,1)-Vr(j,1)).^2+(Vnew(:,2)-Vr(j,2)).^2);
        a=find(Dist<r);
        if ~isempty(a)
            for k=1:length(a)
                id=a(k);
                m=m+1;
                FI(m,1)=Vnew(id,3);
                FI(m,2)=Vr(j,3);
                Vr(j,5)=1;
            end
        end
    end
    FI(:,3)=FI(:,2)./FI(:,1);
    IX=Vr(:,5)==0 & Vr(:,1)>0;
    Fr=Vr(IX,3);
    str=strcat(pathname,filebase,'-',int2str(i),'_result.mat');
    save(str,'FI','Fr');
    str=strcat(pathname,filebase,'-',int2str(i),'_result.xlsx');
    if exist(str,'file')
        delete(str);
    end
    xlswrite(str,FI,['A1:C',num2str(m)]);
    xlswrite(str,Fr,['D1:D',num2str(length(Fr))]);
end