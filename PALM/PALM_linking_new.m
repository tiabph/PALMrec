function VH=PALM_linking_new(VH)

neighbor=VH.parameter.linking.neighborhood;
gap=VH.parameter.linking.gap;

set(VH.text10,'string','Linking particles ...')
mywaitbar(0,VH.axes2,'');
molnumber=0;
VesicleInfo=cell(VH.ImageNumber,1);
n=length(VH.V);        
for i=1:n  
    VesicleInfo{i}=VH.V{i};
    VesicleInfo{i}(:,4)=0;
    molnumber=molnumber+length(VesicleInfo{i}(:,4));
    if floor(i/n*100)>floor((i-1)/n*100) 
        mywaitbar(i/n,VH.axes2,[num2str(floor(i/n*100)),'%']);
    end
end      

nn=0;
tracksFinal=cell(molnumber,1);
mywaitbar(0,VH.axes2,'');
for j=1:VH.ImageNumber
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
            stopframe=min(j+gap,VH.ImageNumber);
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
                        stopframe=min(stopframe+gap,VH.ImageNumber);
                    end                    
                end   
                k=k+1;
            end
        end
    end
    if floor(j/VH.ImageNumber*100)>floor((j-1)/VH.ImageNumber*100) 
    mywaitbar(j/VH.ImageNumber,VH.axes2,[num2str(floor(j/VH.ImageNumber*100)),'%']);
    end
end   
        
trackInfo=cell(nn,1);
mywaitbar(0,VH.axes2,'');
n=length(tracksFinal);
for i=1:n
    if ~isempty(tracksFinal{i})
        trackInfo{i}(:,1)=tracksFinal{i}(:,5);
        trackInfo{i}(:,2:4)=tracksFinal{i}(:,1:3);
    else
        mywaitbar(i/i,VH.axes2,[num2str(floor(i/i*100)),'%']);
        break;
    end
    if floor(i/n*100)>floor((i-1)/n*100) 
        mywaitbar(i/n,VH.axes2,[num2str(floor(i/n*100)),'%']);
    end
end
DV=cell2struct(trackInfo,'trackInfo',2);
VH.DV=DV;