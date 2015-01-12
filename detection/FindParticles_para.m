function V=FindParticles_para(parameter,data_w,v1,v2,v3,w1,w2,A)
parameter.detection.size=0;
if parameter.detection.size>0
    BW=logical(data_w>0);
    [L,num]=bwlabel(BW);
    STATS=regionprops(L,'Area');
    for i=1:num
        if STATS(i).Area>parameter.detection.size
            L(L==i)=1;
        else
            L(L==i)=0;
        end
    end
    data_w=data_w.*L;
end
V=[];
if v1
    if parameter.detection.Cell_area==1
        [W2 W3 A3]=waveletTransform(A,1,3);
        threshold=parameter.detection.boundary_threshold*mean2(A3);
        level=threshold/65535;   
        BW_area=AreaDetect(uint16(A3),level);   
        parameter.detection.Cell_area=BW_area;
    end
    BW_area=parameter.detection.Cell_area;
    data1=data_w.*BW_area;
    if v2
        w=w1;
        ImgMax=locmax2d(data1,[w w]);
        [V(:,2) V(:,1) V(:,3)]=find(ImgMax);
    elseif v3
        w=w2;
        V=weightedcentrid(data1,w);
    end
else   
    ROI=parameter.detection.ROI;
    if ~isempty(ROI)
        data1=data_w.*ROI;
    else
        data1=data_w;
    end
    if v2
        w=w1;
        ImgMax=locmax2d(data1,[w w]);
        [V(:,2) V(:,1) V(:,3)]=find(ImgMax);
    elseif v3
        w=w2;
        V=weightedcentrid_modified(data1,w);
    end
end