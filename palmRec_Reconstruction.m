function databuf = palmRec_Reconstruction(databuf, param)
    recType = param.reconstruction.recType;
    amp = param.reconstruction.amp;
    fitl = param.fitting.fitl;
    
    xmax = (databuf.width +1) * amp;
    ymax = (databuf.height +1) * amp;
    databuf.recimg = zeros(ymax, xmax, 'uint16');
    if(isfield(databuf, 'linkInfo'))
        fitInfo = databuf.linkInfo;
        fitCnt = databuf.linkCnt;
    elseif(isfield(databuf, 'fitInfo'))
        fitInfo = databuf.fitInfo;
        fitCnt = databuf.fitCnt;
    else
        warning('no fit data!');
        return
    end
    
    if(strcmp(recType,'DOT'))
        for i = 1:fitCnt
            cx = round(fitInfo(i,1) * amp) +1;
            cy = round(fitInfo(i,2) * amp) +1;
            if(cx>0 && cy>0 && cx<=xmax && cy<=ymax ...
                    && fitInfo(i,4)<2.5 && fitInfo(i,4)~=1.5 )
                databuf.recimg(cy,cx) = databuf.recimg(cy,cx) + 1;  
            end
        end
    else
        warning('Unrecongnized recType, need "DOT"');
    end
end