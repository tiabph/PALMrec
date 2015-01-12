function V=FindParticles(handles,data_w)
handles.parameter.detection.size=0;
if handles.parameter.detection.size>0
    BW=logical(data_w>0);
    [L,num]=bwlabel(BW);
    STATS=regionprops(L,'Area');
    for i=1:num
        if STATS(i).Area>handles.parameter.detection.size
            L(L==i)=1;
        else
            L(L==i)=0;
        end
    end
    data_w=data_w.*L;
end
V=[];
if get(handles.checkbox1,'value')   
    if handles.parameter.detection.Cell_area==1
        [W2 W3 A3]=waveletTransform(handles.A(:,:,1),1,3);
        threshold=handles.parameter.detection.boundary_threshold*mean2(A3);
        level=threshold/65535;   
    %     BW_area=im2bw(uint16(A3),level); 
    %     BW_area=bwareaopen(BW_area,1000);
        BW_area=AreaDetect(uint16(A3),level);   
        B=bwboundaries(BW_area,'noholes');
        handles.boundary=B{1};
        handles.parameter.detection.Cell_area=BW_area;
    end
    BW_area=handles.parameter.detection.Cell_area;
    data1=data_w.*BW_area;
    if get(handles.radiobutton3,'value')
        w=str2double(get(handles.edit8,'string'));
        ImgMax=locmax2d(data1,[w w]);
        [V(:,2) V(:,1) V(:,3)]=find(ImgMax);
    elseif get(handles.radiobutton4,'value')
        w=str2double(get(handles.edit9,'string'));
        V=weightedcentrid(data1,w);
    end
else   
    handles.ROI=handles.parameter.detection.ROI;
    if ~isempty(handles.ROI)
        data1=data_w.*handles.ROI;
    else
        data1=data_w;
    end
    if get(handles.radiobutton3,'value')
        w=str2double(get(handles.edit8,'string'));
        ImgMax=locmax2d(data1,[w w]);
        [V(:,2) V(:,1) V(:,3)]=find(ImgMax);
    elseif get(handles.radiobutton4,'value')
        w=str2double(get(handles.edit9,'string'));
        V=weightedcentrid(data1,w);
    end
end