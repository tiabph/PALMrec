function databuf = palmRec_CreateSubImages(databuf, param)
    fitl = param.fitting.fitl;
    
    h=waitbar(0, 'Generating sub-image');
    databuf.ROIfit=zeros(2*fitl+1, 2*fitl+1, databuf.pointCnt, 'single');
    for i=1:databuf.pointCnt
        n=databuf.detectionBuf(i,1);
        cx=round(databuf.detectionBuf(i,2));
        cy=round(databuf.detectionBuf(i,3));
        up=cy-fitl; 
        bottom=cy+fitl; 
        left=cx-fitl;
        right=cx+fitl;   
        if up>=1 && left>=1 && bottom<=databuf.height  && right<=databuf.width  && n<databuf.imglen             
            frame=n; 
            I=double(databuf.img(up:bottom,left:right,frame));
            I=I-min(I(:));
            databuf.ROIfit(:,:,i)=single(I);            
        end

        if (floor(i/databuf.pointCnt*100))>(floor((i-1)/databuf.pointCnt*100))
            waitbar(i/databuf.pointCnt, h);
        end
    end
    close(h);
end