a=length(handles.ROI);
[row column]=size(handles.A(:,:,1));
% row=128;
% column=128;
fitl=3;
B=zeros(row,column);
C=zeros(row,column);
D=zeros(row,column);
% [X Y]= meshgrid(1:7,1:7);
% ROI=1000*(exp(-0.5*(X-4).^2./0.3-0.5*(Y-4).^2./0.3));
flag=1;
for k=1:200
    k
nn=0;
for j=flag:a
% m=length(VH.tracksFinalNew{j});
if ~isempty(handles.fitInfo{j})
st=handles.trackInfo{j}(1,1);
x=round(handles.fitInfo{j}(1,1));
y=round(handles.fitInfo{j}(1,2));
sx=handles.fitInfo{j}(1,4);
pN=handles.fitInfo{j}(1,7);
deltaR=handles.fitInfo{j}(1,6);
if st==k && x>fitl && y>fitl && x<=row-fitl && y<=column-fitl %&& sx<2 %&& sx>0.8 && pN>300 && deltaR<50
B(x-fitl:x+fitl,y-fitl:y+fitl)=B(x-fitl:x+fitl,y-fitl:y+fitl)+handles.ROI{j};
% C(x-4:x+4,y-4:y+4)=C(x-4:x+4,y-4:y+4)+handles.ROIfit{j};
C(x,y)=C(x,y)+handles.ROI{j}(fitl+1,fitl+1);
D(x,y)=D(x,y)+1000;
handles.B(:,:,k)=B;
handles.C(:,:,k)=C;
handles.D(:,:,k)=D;
nn=nn+1;
flag=j;
end
end
end
if nn==0
    handles.B(:,:,k)=handles.B(:,:,k-1);
    handles.C(:,:,k)=handles.C(:,:,k-1);
    handles.D(:,:,k)=handles.D(:,:,k-1);
end
end
tiffwrite(handles.B,'F:\B.tif');
tiffwrite(handles.C,'F:\C.tif');
tiffwrite(handles.D,'F:\D.tif');

