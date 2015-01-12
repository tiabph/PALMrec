function databuf = palmRec_DriftCorrection(databuf, param)
    type = param.drift.type;
    if strcmp(type,'file')
        path = param.drift.path;
        file = param.drift.file;
        I = tiffread([path file]);
        framenum = I(3);
        s=size(I);
        fitl = (s(1)-1)/2;
        [X Y]= meshgrid(1:s(1),1:s(2));
        xdata=[X,Y];
        B=zeros(framenum,2);
        h = waitbar(0,'Please wait...');
        for i=1:framenum
            Fd=double(I(:,:,i));
            backg=min(Fd(:));
            peak=max(Fd(:));
            [cy cx]=find(Fd==peak);
            cy=mean(cy);
            cx=mean(cx);
            xp=[cx,cy,1.5,1.5,backg,peak];
            options = optimset('Display','off','MaxFunEvals',1e7,'MaxIter',1000,'TolFun',0.01,'Algorithm','levenberg-marquardt');
            [lp,resnorm,residual,exitflag]=lsqcurvefit(@GaussianKernel,xp,xdata,Fd,[],[],options);
            B(i,2)=lp(1);
            B(i,1)=lp(2);
            [yy xx]=weight_centrid(Fd,s(1));
            if abs(yy-cy)>1 || abs(xx-cx)>1
                yy=cy;
                xx=cx;
            end
            if abs(B(i,1)-yy)>2 || abs(B(i,2)-xx)>2
                B(i,1)=yy;
                B(i,2)=xx;
            end
            if (floor(i/framenum*100))>(floor((i-1)/framenum*100))
                waitbar(i/framenum);
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
        B(:,1)=smooth(B(:,1),50,'rlowess');
        B(:,2)=smooth(B(:,2),50,'rlowess');
        B(:,1)=B(:,1)-B(1,1);
        B(:,2)=B(:,2)-B(1,2);
        databuf.driftInfo = B(:,[2 1]);
    end
end

function f=GaussianKernel(xp,xdata)
    [a b]=size(xdata);
    Xd=xdata(:,1:b/2);
    Yd=xdata(:,b/2+1:b);
    f=xp(6)*(exp(-0.5*(Xd-xp(1)).^2./(xp(3)^2)-0.5*(Yd-xp(2)).^2./(xp(4)^2)))+xp(5);
end


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
end