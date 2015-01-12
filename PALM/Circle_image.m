A=zeros(64,64);
for i=1:64
    for j=1:64
        dist=sqrt((i-32)^2+(j-32)^2);
        if round(dist)==8  || round(dist)==16 || round(dist)==24 %|| round(dist)
%         if (dist<=4 && dist>=4) || (dist<=6 && dist>=6) || (dist<=8 && dist>=8) || (dist<=10 && dist>=10) || (dist<=12 && dist>=12) 
            %|| (dist<=30 && dist>=28) || (dist<=34 && dist>=32) ...
             %   || (dist<=38 && dist>=36) || (dist<=42 && dist>=40)  || (dist<=46 && dist>=44) || (dist<=50 && dist>=48) 
            A(i,j)=1;
        end
    end
end
% figure;imshow(A,[])
tiffwrite(A*255)
