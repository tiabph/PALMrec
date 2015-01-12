clear
clc
peak=1000;
fitl=4;
[Xd,Yd] = meshgrid(-fitl:fitl,-fitl:fitl);
xdata=[Xd,Yd];
sx=1.2;
cx=0.25;
cy=0;
ROI1=peak*(exp(-0.5*(Xd-cx).^2./(sx^2)-0.5*(Yd-cy).^2./(sx^2)));
cx=0.5;
cy=0;
ROI2=peak*(exp(-0.5*(Xd-cx).^2./(sx^2)-0.5*(Yd-cy).^2./(sx^2)));
cx=0.75;
cy=0;
ROI3=peak*(exp(-0.5*(Xd-cx).^2./(sx^2)-0.5*(Yd-cy).^2./(sx^2)));
I=zeros(64,32,1000);
Y=5:59;
x=16;
for i=1:1000
    i
    N=round(54*rand(1)+1);
    y1=Y(N);
    I(y1-4:y1+4,x-4:x+4,i)=I(y1-4:y1+4,x-4:x+4,i)+ROI1;
    N=round(54*rand(1)+1);
    y2=Y(N);
    while abs(y2-y1)<9
        N=round(54*rand(1)+1);
        y2=Y(N);
    end
    I(y2-4:y2+4,x-4:x+4,i)=I(y2-4:y2+4,x-4:x+4,i)+ROI2;
    N=round(54*rand(1)+1);
    y3=Y(N);
    while abs(y3-y2)<9 || abs(y3-y1)<9
        N=round(54*rand(1)+1);
        y3=Y(N);
    end
    I(y3-4:y3+4,x-4:x+4,i)=I(y3-4:y3+4,x-4:x+4,i)+ROI3;
end
tiffwrite(I,'F:\12.tif');
