addpath('..\util_cpu\');
drift_gap = 1;
drift_dist = 0.5;
[re sp] = LinkPoints(databuf.fitInfo, [1 2 9], size(databuf.img),[drift_gap drift_dist]);

driftbuf = zeros(databuf.imglen, 3);
for m=1:size(sp,1)
    if(sp(m,2) >1) %linked points
        plist = zeros(1, sp(m,2));
        cp = sp(m,1);
        pcnt = 1;
        plist(1) = cp-1;
        while(1)
            np = re(cp);
            if(np<0)
                break;
            end
            pcnt = pcnt +1;
            plist(pcnt) = np+1;
            cp=np;
        end
        
        if(pcnt ~= sp(m,2))
            disp('link err');
        end
        
        for n=1:pcnt-1
            dx = databuf.fitInfo(plist(n+1), 1) - databuf.fitInfo(plist(n), 1);
            dy = databuf.fitInfo(plist(n+1), 2) - databuf.fitInfo(plist(n), 2);
            fs = databuf.fitInfo(plist(n),9);
            fe = databuf.fitInfo(plist(n+1),9);
            df = fe-fs;
            dx = dx/df;
            dy = dy/df;
            for p = (fs+1):fe
                driftbuf(p,1) = driftbuf(p,1) + dx;
                driftbuf(p,2) = driftbuf(p,2) + dy;
                driftbuf(p,3) = driftbuf(p,3) + 1;
            end
        end
    end
end

drift_d = driftbuf(:, 1:2);
drift_d(2:end,1) = drift_d(2:end,1) ./ driftbuf(2:end,3);
drift_d(2:end,2) = drift_d(2:end,2) ./ driftbuf(2:end,3);

drift_i = drift_d;
for m=2:size(drift_i,1)
    drift_i(m, 1) = drift_i(m, 1) + drift_i(m-1, 1);
    drift_i(m, 2) = drift_i(m, 2) + drift_i(m-1, 2);
end

