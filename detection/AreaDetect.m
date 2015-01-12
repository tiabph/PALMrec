%set the largest area above threshold "level(0-1)" to 1, return the binary
%image BW_Z
function BW_Z = AreaDetect(AA,level)

[row,column] = size(AA);
BW_A = im2bw(AA,level);
[L,num]=bwlabel(BW_A);
maxnum=1;
[r,c]=find(L==1);
maxarea=size(r,1);
BW_Z=BW_A;
if(num>1)
    for i=2:num
        [r,c]=find(L==i);
        if(size(r,1)>maxarea)
            maxarea=size(r,1);
            maxnum=i;
        end
    end   
    [r,c]=find(L==maxnum);
    BW_Z=zeros(row,column);
    for k=1:maxarea
    BW_Z(r(k),c(k))=1;
    end
end
