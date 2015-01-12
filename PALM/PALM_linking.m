function VH=PALM_linking(VH,neighbor,gap)

    for i=1:length(VH.V)        
        VesicleInfo{i}=VH.V{i};
        VesicleInfo{i}(:,4)=0;
    end    
    nn=0;
    set(VH.text10,'string','Linking particles ...')
    mywaitbar(0,VH.axes2,'');
    for j=1:VH.ImageNumber-1       
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
            for k=j+1:VH.ImageNumber
                TNx=VesicleInfo{k}(:,1);
                TNy=VesicleInfo{k}(:,2);
                Dist=(TNx-xx).^2+(TNy-yy).^2;
                if neighbor==8
                    FI=find(Dist<3,1);
                elseif neighbor==4
                    FI=find(Dist<2,1);
                end                
                if FI>0 
                    if VesicleInfo{k}(FI,4)==0
                    VesicleInfo{k}(FI,4)=1;
                    tracksFinal{nn}(ii,1:4)=VesicleInfo{k}(FI,1:4);
                    xx=tracksFinal{nn}(ii,1);
                    yy=tracksFinal{nn}(ii,2);
                    tracksFinal{nn}(ii,5)=k;
                    ii=ii+1;
                    end
                end
            end
            end
        end
        mywaitbar(j/(VH.ImageNumber-1),VH.axes2,[num2str(floor(j/(VH.ImageNumber-1)*100)),'%']);
    end   
    
    tracksNew=tracksFinal;
    tracksFinal=[];
    m=length(tracksNew);
    nn=1;
    set(VH.text10,'string','Linking tracks ...')
    mywaitbar(0,VH.axes2,'');
    for i=1:m
        l=length(tracksNew{i}(:,1));
        if l==1
           tracksFinal{nn}(:,1:5)=tracksNew{i}(:,1:5);       % only one frame then one track
           nn=nn+1;
        end
        if l>1
            sp=[];
            ii=1;
            for j=2:l
                if (tracksNew{i}(j,5)-tracksNew{i}(j-1,5)>gap)
                    sp(ii,1)=j;                              % find break point,ie. break for more than 2 frames
                    ii=ii+1;
                end
            end
            n=length(sp);
            if n==0  % no gap
                 tracksFinal{nn}(:,1:5)=tracksNew{i}(:,1:5);     
                 nn=nn+1;
            end
            if n>0  
                tracksFinal{nn}(:,1:5)=tracksNew{i}(1:sp(1,1)-1,1:5); 
                nn=nn+1;  
                if n>1
                for k=1:n-1
                    tracksFinal{nn}(:,1:5)=tracksNew{i}(sp(k,1):sp(k+1,1)-1,1:5); 
                    nn=nn+1;            
                end  
                end
                tracksFinal{nn}(:,1:5)=tracksNew{i}(sp(n,1):l,1:5); 
                nn=nn+1;  
            end
        end
        mywaitbar(i/m,VH.axes2,[num2str(floor(i/m*100)),'%']);
    end
        
    DV=[];
    for i=1:length(tracksFinal)
        DV(i,1).trackInfo(:,1)=tracksFinal{i}(:,5);
        DV(i,1).trackInfo(:,2:4)=tracksFinal{i}(:,1:3);
    end
    VH.DV=DV;