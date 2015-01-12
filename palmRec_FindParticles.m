function databuf = palmRec_FindParticles(databuf, param)
    threshold = param.detection.threshold;
    type = param.detection.type;
    windowWidth = param.detection.windowWidth;
    
    databuf.detectionResult = cell(databuf.imglen, 1);
    h = waitbar(0,'Detection progress');
    pointCnt = 0;
    for m=1:databuf.imglen
        data_w =Detection(databuf.img(:,:,m), threshold, type);
        databuf.detectionResult{m} = weightedcentrid(data_w,windowWidth);
        pointCnt = pointCnt + size(databuf.detectionResult{m},1);

        if (floor(m/databuf.imglen*100))>(floor((m-1)/databuf.imglen*100))
            waitbar(m/databuf.imglen, h);
        end
    end
    databuf.pointCnt = pointCnt;
    close(h);
    databuf.detectionBuf = zeros(pointCnt, 4);
    ts=1;
    for m=1:databuf.imglen
        temp = databuf.detectionResult{m};
        tempLen = size(temp,1);
        temp = cat(2, ones(tempLen,1).*m, temp);
        databuf.detectionBuf(ts:ts+tempLen-1, :) = temp;
        ts = ts+tempLen;
    end
end