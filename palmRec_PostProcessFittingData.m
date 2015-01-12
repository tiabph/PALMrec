function databuf = palmRec_PostProcessFittingData(databuf, param)
    h=waitbar(0, 'post-process fit info');
    fitl = param.fitting.fitl;
    factor = param.fitting.factor;
    gain = param.fitting.gain;
    a=param.fitting.pixelsize;
    
    fitInfo = databuf.fitInfo_raw;
    databuf.fitInfo=zeros(databuf.pointCnt,9);

    fitcnt = 0;
    bkgBuf = zeros(databuf.pointCnt,2,(fitl*2+1)^2);
    tempbuf = zeros(databuf.pointCnt,2);
    for i=1:databuf.pointCnt
        cx=round(databuf.detectionBuf(i,2));
        cy=round(databuf.detectionBuf(i,3));
        rcx=round(cx);
        rcy=round(cy);
        n=databuf.detectionBuf(i,1);
        up=round(cy)-fitl; 
        bottom=round(cy)+fitl; 
        left=round(cx)-fitl;
        right=round(cx)+fitl;       
        if up>=1 && left>=1 && bottom<=databuf.height  && right<=databuf.width  && n<databuf.imglen ...
            && fitInfo(i,1)~= fitl+1 && fitInfo(i,2)~= fitl+1 ...
            && fitInfo(i,1)~= 0 && fitInfo(i,2)~= 0 ...
            && fitInfo(i,1)~= fitl*2+1 && fitInfo(i,2)~= fitl*2+1
    % ---------- position ----------
            cx=rcx+fitInfo(i,1)-fitl-1;
            cy=rcy+fitInfo(i,2)-fitl-1; 
    % ---------- peak ----------
            peak=fitInfo(i,3);
    % ---------- sigma x and y ----------
            sx=fitInfo(i,5);
            sy=fitInfo(i,5);
    % ---------- photon number ----------
            Nm=fitInfo(i,3);
    % ---------- frame number ----------
            first_frame = n;
            if(first_frame >1)
                I1=(databuf.img(rcy-fitl:rcy+fitl,rcx-fitl:rcx+fitl,first_frame-1)); 
            else
                I1 = NaN;
            end
            if(first_frame < databuf.pointCnt)
                I2=(databuf.img(rcy-fitl:rcy+fitl,rcx-fitl:rcx+fitl,first_frame+1));
            else
                I2 = NaN;
            end
            
%             sd = min(std(I1(:)),std(I2(:)));
%             b=sd*factor/gain;  
            s=sx*a;    
            
            
%             deltaR=sqrt((s^2+a^2/12)/Nm+8*pi*s^4*b^2/a^2/Nm^2);
            fitcnt = fitcnt+1;
            bkgBuf(fitcnt,1,:) = I1(:);
            bkgBuf(fitcnt,2,:) = I2(:);
            tempbuf(fitcnt,1) = (s^2+a^2/12)/Nm;
            tempbuf(fitcnt,2) = 8*pi*s^4/a^2/Nm^2;
            databuf.fitInfo(fitcnt,:)=[cx cy peak sx sy 0 Nm 0 first_frame];
        end
        
        if (floor(i/databuf.pointCnt*100))>(floor((i-1)/databuf.pointCnt*100))
            waitbar(i/databuf.pointCnt, h);
        end
    end

    databuf.fitInfo = databuf.fitInfo(1:fitcnt,:);
    bkgBuf = bkgBuf(1:fitcnt,:,:);
    tempbuf = tempbuf(1:fitcnt,:);
    b = min(std(bkgBuf,0,3),[],2) *factor/gain;
    deltaR=sqrt(tempbuf(:,1)+tempbuf(:,2).*b.^2);
    databuf.fitInfo(:,6) = deltaR;
    databuf.fitInfo(:,8) = b;
    
    databuf.fitCnt = fitcnt;
    close(h);
end