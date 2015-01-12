% clear
clc
pathname='G:\zms\20130915_meos3.2-lifeact\4\';
filebase='cell4';
N=4;    % stack number
r=2;    % distance (pixel)
t=0;    % time
for i=1:1
    i
    FI=[];
    m=0;
    str=strcat(pathname,filebase,'_l_',int2str(i),'_Reconstruction.mat');
    P=importdata(str);
    Vl=P.V;
    str=strcat(pathname,filebase,'_r_',int2str(i),'_Reconstruction.mat');
    P=importdata(str);
    Vr=P.V;
    Vr(:,5)=0;
    P=[];
    nl=length(Vl(:,1));
    nr=length(Vr(:,1));
    for j=1:nr
        IX=abs(Vl(:,3)-Vr(j,3))<=t;
        Vnew=Vl(IX,:);
        Dist=sqrt((Vnew(:,1)-Vr(j,1)).^2+(Vnew(:,2)-Vr(j,2)).^2);
        a=find(Dist<r);
        if ~isempty(a)
            for k=1:length(a)
                id=a(k);
                m=m+1;
                FI(m,1)=Vnew(id,4);
                FI(m,2)=Vr(j,4);
                Vr(j,5)=1;
            end
        end
    end
    FI(:,3)=FI(:,2)./FI(:,1);
    IX=Vr(:,5)==0;
    Fr=Vr(IX,4);
    str=strcat(pathname,filebase,'-',int2str(i),'_result.xlsx');
    if exist(str,'file')
        delete(str);
    end
    xlswrite(str,FI,['A1:C',num2str(m)]);
    xlswrite(str,Fr,['D1:D',num2str(length(Fr))]);
end