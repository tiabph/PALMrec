% clear
% clc
%     
%     I=tiffread('E:\test.tif');
%     I=single(I);
%         kernel{1}=[1/16,1/4,3/8,1/4,1/16];
% %     kernel{2}=[1/16,0,1/4,0,3/8,0,1/4,0,1/16];
% %     kernel{3}=[1/16,0,0,0,1/4,0,0,0,3/8,0,0,0,1/4,0,0,0,1/16];
%     kernel{2}=[1/16,1/4,3/8,1/4,1/16];
%     kernel{3}=[1/16,1/4,3/8,1/4,1/16];
%     
%     gfor k=1:10
%     A=I(:,:,k);
%     [row column]=size(A);
%     A0=A;
%     A=gzeros(row,column,10);
%     W=gzeros(row,column,10);
% 
%     for i=1:3
%         if i==1
%             I=A0;
%             kernel1=kernel{1};
%             [row column]=size(I);
%             extend=(length(kernel1)-1)/2;
% %             I=padarray(I,[extend extend],'symmetric');
%             I1=conv2(kernel1,kernel1',I,'same');
% %             A(:,:,1)=I1(extend+1:extend+row,extend+1:extend+column);  
%             A(:,:,1)=I1;
%             W(:,:,i)=A0-A(:,:,1);
%         else
%             kernel1=kernel{i};
%             I=A(:,:,i-1);
%             [row column]=size(I);
%             extend=(length(kernel1)-1)/2;
% %             I=padarray(I,[extend extend],'symmetric');
%             I1=conv2(kernel1,kernel1',I,'same');
% %             A(:,:,i)=I1(extend+1:extend+row,extend+1:extend+column);
%             A(:,:,i)=I1;
%             W(:,:,i)=A(:,:,i-1)-A(:,:,i);
%         end
%     end    
%     
%     for i=1:3
%         deta=median(median(abs(W{i}-median(median(W{i})).*ones(row,column))));
%         t=coef*deta/0.67;
%         W(:,:,i)=(W(:,:,i)>=t).*W(:,:,i);
%     end
%           
% 
%     varargout{1}=W(:,:,1);
%     varargout{2}=W(:,:,2);
%     varargout{3}=W(:,:,3);
%     varargout{4}=A(:,:,3);
%     
%     gend
flag=3;
for i=1:flag
    flag=flag+1;
    flag=min(flag,10);
    i
end
