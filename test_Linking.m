%% test script for Linking

param.linking.gap = 1;
param.linking.neighbor = 8;

%matlab linking
databuf = palmRec_Linking(databuf, param);

%cpu linking
databuf = palmRec_Linking_cpu(databuf, param);

%% check link result
linkBuf_raw = databuf.linkResult;
linkBuf_re = databuf.linkData_chain;
linkBuf_sp = databuf.linkData_sp;

linkBuf_cpu = cell(size(linkBuf_sp,1),1);

for m=1:size(linkBuf_sp,1)
    cp = linkBuf_sp(m,1)+1;
    linklist = zeros(linkBuf_sp(m,2),1);
    linklist(1) = cp;
    cnt = 1;
    while(linkBuf_re(cp) >=0)
        np = linkBuf_re(cp) +1;
        cnt=cnt+1;
        linklist(cnt) = np;
        cp=np;
    end
    linkBuf_cpu{m} = linklist;
end

%% check link flag
linkFlag_raw = zeros(databuf.fitCnt, 1);
linkFlag_cpu = zeros(databuf.fitCnt, 1);

for m=1:length(linkBuf_raw)
    linklist = linkBuf_raw{m};
    for n=1:(length(linklist)-1)
        linkFlag_raw(linklist(n)) = linklist(n+1);
    end
end

for m=1:length(linkBuf_cpu)
    linklist = linkBuf_cpu{m};
    for n=1:(length(linklist)-1)
        linkFlag_cpu(linklist(n)) = linklist(n+1);
    end
end

cmpresult.same = sum(linkFlag_cpu == linkFlag_raw);
cmpresult.cpu = sum( (linkFlag_cpu ~= linkFlag_raw) & (linkFlag_cpu ==0) );
cmpresult.raw = sum( (linkFlag_cpu ~= linkFlag_raw) & (linkFlag_raw ==0) );
cmpresult.diff = sum( linkFlag_cpu ~= linkFlag_raw ) - cmpresult.cpu - cmpresult.raw;

%% output report
disp('---------- test Linking report ----------');
fprintf(1,'total points: %d\nsame : %d (%.2f%%)\ncpu : %d (%.2f%%)\nraw : %d (%.2f%%)\ndiff : %d (%.2f%%)\n', ...
                    databuf.fitCnt, ...
                    cmpresult.same, cmpresult.same/databuf.fitCnt *100, ...
                    cmpresult.cpu, cmpresult.cpu/databuf.fitCnt *100, ...
                    cmpresult.raw, cmpresult.raw/databuf.fitCnt *100, ...
                    cmpresult.diff, cmpresult.diff/databuf.fitCnt *100 ...
                    );
disp('------------- end of report -------------');        
