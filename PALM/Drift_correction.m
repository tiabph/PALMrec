clear
clc
sigma = 1.5;
pathname = 'F:\data\20141022_orai1\gt\cell2\1\';
filebase='cell1_beads';
str = strcat(pathname,filebase,'.tif');
[I, n] = tiffread(str);
wid=size(I,1);
r = size(I,1);
fitl = (r-1)/2;
[X Y]= meshgrid(1:wid,1:wid);
xdata=[X,Y];
B=zeros(n,2);
h = waitbar(0,'Please wait...');
for i=1:n
    Fd=double(I(:,:,i));
    backg=min(Fd(:));
    peak=max(Fd(:));
    [cy cx]=find(Fd==peak);
    cy=mean(cy);
    cx=mean(cx);
    xp=[cx,cy,1.5,1.5,backg,peak];
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
        waitbar(i/n);
    end     
end
% for i=2:length(B)-1
%     if abs(B(i,1)-B(i-1,1))>3 && abs(B(i-1,1)-B(i+1,1))<0.5
%         B(i,1)=(B(i-1,1)+B(i+1,1))/2;
%     end
%     if abs(B(i,2)-B(i-1,2))>3 && abs(B(i-1,2)-B(i+1,2))<0.5
%         B(i,2)=(B(i-1,2)+B(i+1,2))/2;
%     end
% end
close(h);
pause(0.1);
B(:,1)=smooth(B(:,1),50,'rlowess');
B(:,2)=smooth(B(:,2),50,'rlowess');
B(:,1)=B(:,1)-B(1,1);
B(:,2)=B(:,2)-B(1,2);
figure;plot(B(:,1),'b');hold on;plot(B(:,2),'r');
str=strcat(pathname,filebase,'.mat');
save(str,'B');