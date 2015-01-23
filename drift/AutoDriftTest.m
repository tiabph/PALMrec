addpath('..\util_cpu\');
drift_gap = 50;
drift_dist = 0.5;
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
        weight = 1;
        
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

% fitInfo = databuf.fitInfo_raw;
fitInfo(:,1) = fitInfo(:,1) - drift_i(fitInfo(:,9),1);
fitInfo(:,2) = fitInfo(:,2) - drift_i(fitInfo(:,9),2);

%% display
figure(1)
plot(drift_i(:,1),'r');
hold on
plot(drift_i(:,2),'g');
hold off
legend('x','y');

% figure(2)
% plot(databuf.fitInfo(:,1),databuf.fitInfo(:,2),'b.')
% hold on
% plot(fitInfo(:,1),fitInfo(:,2),'r.')
% hold off