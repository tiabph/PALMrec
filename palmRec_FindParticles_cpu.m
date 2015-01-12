function databuf = palmRec_FindParticles_cpu(databuf, param)
    threshold = param.detection.threshold;
    dettype = param.detection.type;
    windowWidth = param.detection.windowWidth;
    
    databuf.detectionResult = DWTParticalDetection(databuf.img, threshold, dettype, windowWidth);
    pointCnt = 0;
    for m=1:databuf.imglen
        pointCnt = pointCnt + size(databuf.detectionResult{m},1);
    end
    databuf.pointCnt = pointCnt;

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