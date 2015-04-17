function databuf = palmRec_Filter(databuf, param)
    if ~isfield(param, 'filter')
        return
    end
    maxPositionError = param.filter.max_position_err;
    maxPhoton = param.filter.max_photon_num;
    minPhoton = param.filter.min_photon_num;
    maxSigma = param.filter.max_sigma;
    minSigma = param.filter.min_sigma;
    maxBackground = param.filter.max_backgound;
    maxLife = param.filter.max_lifetime;
    
    info = databuf.linkInfo;
    mask = true(databuf.linkCnt,1);
    %position error
    mask = mask & (info(:,6)<maxPositionError);
    
    %photon number
    mask = mask & (info(:,7)<maxPhoton);
    mask = mask & (info(:,7)>minPhoton);
    
    %sigma
    mask = mask & (info(:,4)<maxSigma);
    mask = mask & (info(:,4)>minSigma);
    mask = mask & (info(:,5)<maxSigma);
    mask = mask & (info(:,5)>minSigma);
    
    %Background
    mask = mask & (info(:,8)<maxBackground);
    
    %life time
    mask = mask & (databuf.linkData_sp(:,2)<maxLife);
    
    databuf.linkInfo = info(mask,:);
    databuf.linkCnt = size(databuf.linkInfo, 1);
end