function [Image,number]=tiffread(FileStr,idx)
%% Reading a tiff image into a row*column*N matrix
% *******************************************************
% [Image,number,ImageNumber]=tiffread();
% [Image,number,ImageNumber]=tiffread(FileStr);
% [Image,number,ImageNumber]=tiffread(FileStr,idx);
% *******************************************************
% Author:ZYD,IBP,CAS 11/30/2009

if nargin==0
    [filename, pathname] = uigetfile( ...
        {'*.tif;*.tiff', 'All TIF-Files (*.tif,*.tiff)'; ...
            '*.*','All Files (*.*)'}, ...
        'Open Image File');
    if isequal([filename,pathname],[0,0])
        return
    else
        FileStr = fullfile(pathname,filename);
    end
end

A=readtif(FileStr,1);
[row,col]=size(A);

ImageInfo = imfinfo(FileStr);

if nargin<2
    h=waitbar(0,'Reading images,please wait...');
%     try
%         info=tifftagsread(FileStr);
%     catch me
        info=imfinfo(FileStr);
%     end
    ImageNumber=length(info);
    Image=uint16(zeros(row,col,ImageNumber));
    for j=1:ImageNumber
        Image(:,:,j)=readtif(FileStr, j);
        if (floor(j/ImageNumber*100))>(floor((j-1)/ImageNumber*100))
            waitbar(j/ImageNumber, h);
        end
    end
    close(h);
    number=ImageNumber;
    pause(eps)
end

if nargin==2
    if length(idx)==1
       Image=uint16(zeros(row,col));
       Image(:,:)=readtif(FileStr, idx); 
       number=1;
    else
       n=idx(2)-idx(1)+1;
       v=idx(1):idx(2);
       Image=uint16(zeros(row,col,n));
       h=waitbar(0,'Reading images,please wait...');
       for j=1:n
           Image(:,:,j)=readtif(FileStr, v(j));
           waitbar(j/n);
       end
       number=n;
       close(h);
       pause(eps)
    end
end
