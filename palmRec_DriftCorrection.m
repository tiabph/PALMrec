function databuf = palmRec_DriftCorrection(databuf, param)
    type = param.drift.type;
    B=[];
    if isempty(param.drift.file)
        databuf.DriftCorrectionFlag = 0;
        return
    end
    
    if strcmp(type,'file')
        path = param.drift.path;
        file = param.drift.file;
        I = TIFFloadFrame_16bit_CPU([path file],[1 -1]);
        %kernel function
%         B = CalDrift_raw(I);    %raw function with gaussian fitting
%         B = CalDrift_weightedcentroid(I);  %weighted centroid 
        B = CalDrift_GPU(I);    %gaussian fitting with gpu
        B=double(B);
        
        B(:,1) = B(:,1) - B(1,1);
        B(:,2) = B(:,2) - B(1,2);
        
%         B(:,1) = smooth(B(:,1), 50, 'rlowess');
%         B(:,2) = smooth(B(:,2), 50, 'rlowess');
        
        B(:,1) = medfilt1(B(:,1), 50);
        B(:,2) = medfilt1(B(:,2), 50);
        
        B(:,1) = B(:,1) - B(1,1);
        B(:,2) = B(:,2) - B(1,2);
        databuf.driftInfo = B;
    end
    
    if( ~isempty(B))%start correction
        databuf.fitInfo_raw = databuf.fitInfo;
        framelist = databuf.fitInfo(:,9);
        x_corr = B(framelist,1);
        y_corr = B(framelist,2);
        databuf.fitInfo(:,1) = databuf.fitInfo(:,1) - x_corr;
        databuf.fitInfo(:,2) = databuf.fitInfo(:,2) - y_corr;
        
        databuf.DriftCorrectionFlag = 1;
    end
end

