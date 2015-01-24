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
    elseif strcmp(type,'auto')
        drift_gap = param.drift.gap;
        drift_dist = param.drift.dist;
        [re sp] = LinkPoints(databuf.fitInfo, [1 2 9], size(databuf.img),[drift_gap drift_dist]);
        fitInfo = databuf.fitInfo;
        driftbuf = zeros(databuf.imglen, 3);
        for m=1:size(sp,1)
            if(sp(m,2) >2) %linked points
                plist = zeros(1, sp(m,2));
                cp = sp(m,1) +1;
                pcnt = 1;
                plist(1) = cp;
                while(1)
                    np = re(cp) +1;
                    if(np<=0)
                        break;
                    end
                    pcnt = pcnt +1;
                    plist(pcnt) = np;
                    cp=np;
                end

                if(pcnt ~= sp(m,2))
                    disp('link err');
                end

                %---------- weight of point ----------
        %         weight = databuf.fitInfo(plist, 3);
        %         weight = sum(weight(:));
                weight = pcnt^2;

                for n=1:pcnt-1
                    dx = fitInfo(plist(n+1), 1) - fitInfo(plist(n), 1);
                    dy = fitInfo(plist(n+1), 2) - fitInfo(plist(n), 2);
                    fs = fitInfo(plist(n),9);
                    fe = fitInfo(plist(n+1),9);
                    df = fe-fs;
                    dx = dx/df;
                    dy = dy/df;
                    for p = (fs+1):fe
                        driftbuf(p,1) = driftbuf(p,1) + dx*weight;
                        driftbuf(p,2) = driftbuf(p,2) + dy*weight;
                        driftbuf(p,3) = driftbuf(p,3) + weight;
                    end
                end
            end
        end

        drift_d = driftbuf(:, 1:2);
        drift_d(2:end,1) = drift_d(2:end,1) ./ driftbuf(2:end,3);
        drift_d(2:end,2) = drift_d(2:end,2) ./ driftbuf(2:end,3);

        drift_d(2:end,1) = smooth(drift_d(2:end,1), 50, 'rlowess');
        drift_d(2:end,2) = smooth(drift_d(2:end,2), 50, 'rlowess');

        drift_i = drift_d;
        for m=2:size(drift_i,1)
            drift_i(m, 1) = drift_i(m, 1) + drift_i(m-1, 1);
            drift_i(m, 2) = drift_i(m, 2) + drift_i(m-1, 2);
        end
        B = drift_i;
        databuf.driftInfo = drift_i;
    end
    
    if( isfield(databuf, 'driftInfo') )%start correction
        databuf.fitInfo_raw = databuf.fitInfo;
        framelist = databuf.fitInfo(:,9);
        x_corr = B(framelist,1);
        y_corr = B(framelist,2);
        databuf.fitInfo(:,1) = databuf.fitInfo(:,1) - x_corr;
        databuf.fitInfo(:,2) = databuf.fitInfo(:,2) - y_corr;
        
        databuf.DriftCorrectionFlag = 1;
    end
end

