function databuf = palmRec_Fitting(databuf, param)
    diary('fitlog.txt');
    diary on;
    blockSize = 8192*8;
    factor = param.fitting.factor;
    gain = param.fitting.gain;
    sigma = param.fitting.sigma;
    N=floor(databuf.pointCnt/blockSize);
    databuf.fitInfo_raw=[];
    for i=1:N+1
        st=(i-1)*blockSize+1;
        et=i*blockSize;
        et=min(et,databuf.pointCnt);
        ROI=databuf.ROIfit(:,:,st:et)*factor/gain;
%         tic
%         [Y X N BG S CRLBy CRLBx CRLBn CRLBb CRLBs LogL]=mGPUgaussMLE(single(ROI),sigma,20,2);
%         toc
        [ret, Y X N BG S CRLBy CRLBx CRLBn CRLBb CRLBs LogL] = evalc('mGPUgaussMLE(single(ROI),sigma,20,2);');
        databuf.fitInfo_raw(st:et,1)=X+1;
        databuf.fitInfo_raw(st:et,2)=Y+1;
        databuf.fitInfo_raw(st:et,3)=N;
        databuf.fitInfo_raw(st:et,4)=BG;
        databuf.fitInfo_raw(st:et,5)=S;
    end
%     databuf.fitInfo_raw=fitInfo;
    diary off;
end