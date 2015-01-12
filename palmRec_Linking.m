function databuf = palmRec_Linking(databuf, param)
    gap=param.linking.gap;
    neighbor=param.linking.neighbor;
    imglen = databuf.imglen;
    
    molnumber=0;
    VesicleInfo=cell(databuf.imglen,1);
    indexList = 1:size(databuf.fitInfo, 1);
    for i=1:databuf.imglen
        mask = (databuf.fitInfo(:,9)==i);
        VesicleInfo{i}=databuf.fitInfo(mask,[1 2 9 9]);
        if(~isempty(find(mask,1)))
            VesicleInfo{i}(:,3)=indexList(mask);
            VesicleInfo{i}(:,4)=0;
        end
        molnumber=molnumber+length(VesicleInfo{i}(:,4));
    end 
    
    nn=0;
    tracksFinal=cell(molnumber,1);
    for j=1:imglen
        m=length(VesicleInfo{j}(:,1));
        for i=1:m           
            if VesicleInfo{j}(i,4)==0
                nn=nn+1;
                tracksFinal{nn}(1,1:3)=VesicleInfo{j}(i,1:3);
                tracksFinal{nn}(1,4)=1;
                tracksFinal{nn}(1,5)=j;
                ii=2;
                xx=VesicleInfo{j}(i,1);
                yy=VesicleInfo{j}(i,2);
                stopframe=min(j+gap,imglen);
                k=j+1;
                while k<=stopframe                
                    TNx=VesicleInfo{k}(:,1);
                    TNy=VesicleInfo{k}(:,2);
                    Dist=(TNx-xx).^2+(TNy-yy).^2;
                    if neighbor==8
                        FI=find(Dist<1.5,1);
                    elseif neighbor==4
                        FI=find(Dist<1,1);
                    elseif neighbor==1
                        FI=find(Dist<0.5,1);
                    end
                    if FI>0 
                        if VesicleInfo{k}(FI,4)==0
                            VesicleInfo{k}(FI,4)=1;
                            tracksFinal{nn}(ii,1:4)=VesicleInfo{k}(FI,1:4);
                            xx=tracksFinal{nn}(ii,1);
                            yy=tracksFinal{nn}(ii,2);
                            tracksFinal{nn}(ii,5)=k;
                            ii=ii+1;
                            stopframe=min(stopframe+gap,imglen);
                        end                    
                    end   
                    k=k+1;
                end
            end
        end
    end  
    
    trackInfo=cell(nn,1);
    n=nn;
    for i=1:n
        if ~isempty(tracksFinal{i})
%             trackInfo{i}(:,1)=tracksFinal{i}(:,5);
%             trackInfo{i}(:,2:4)=tracksFinal{i}(:,1:3);
            trackInfo{i}=tracksFinal{i}(:,3);
        else
            break;
        end
    end
    databuf.linkResult=trackInfo; 
    databuf.linkCnt = n;
    
    linkInfo = zeros(databuf.linkCnt, 9);
    for i = 1:databuf.linkCnt
        plist = databuf.linkResult{i};
        tempinfo = databuf.fitInfo(plist,:);
        linkInfo(i,1) = mean(tempinfo(:,1));%cx
        linkInfo(i,2) = mean(tempinfo(:,2));%cy
        linkInfo(i,3) = sum(tempinfo(:,3));%peak
        linkInfo(i,4) = mean(tempinfo(:,4));%sx
        linkInfo(i,5) = mean(tempinfo(:,5));%sy
        linkInfo(i,6) = mean(tempinfo(:,6));%deltaR
        linkInfo(i,7) = sum(tempinfo(:,7));%Nm
        linkInfo(i,8) = mean(tempinfo(:,8));%b
        linkInfo(i,9) = mean(tempinfo(1,9));%first frame
    end
    databuf.linkInfo = linkInfo;
    
end