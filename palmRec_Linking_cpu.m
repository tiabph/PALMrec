function databuf = palmRec_Linking_cpu(databuf, param)
    gap=param.linking.gap;
    imglen = databuf.imglen;
    if(isfield(param.linking, 'neighbor'))
        neighbor=param.linking.neighbor;
        if(neighbor ==1)
            dist = 0.5;
        elseif(neighbor == 4)
            dist = 1.0;
        else
            dist = 1.5;
        end
    else
        dist = param.linking.dist;
    end
    
    [re sp] = LinkPoints(databuf.fitInfo, [1 2 9], size(databuf.img),[gap dist]);
    databuf.linkData_chain = re;
    databuf.linkData_sp = sp;
    databuf.linkInfo = PostProcessLinkData(databuf.fitInfo, re, sp);
    databuf.linkCnt = size(databuf.linkInfo ,1);
end