function [W2 W3] = det_DWT_gpu(img)
%wavelet tranform for partical detection

    kernel{1}=[1/16,1/4,3/8,1/4,1/16];
    kernel{2}=[1/16,0,1/4,0,3/8,0,1/4,0,1/16];
    kernel{3}=[1/16,0,0,0,1/4,0,0,0,3/8,0,0,0,1/4,0,0,0,1/16];
    [row column imglen]=size(img);
    threadlen = 1000;
%     A0=single(A);
%     A=cell(1,3);
%     T=cell(1,3);
%     W=cell(1,3);
%     W2=zeros(row, column, imglen);
%     W3=zeros(row, column, imglen);
% 
%     for imgcnt = 1:imglen
%         A0 = double(img(:,:,imgcnt));
%         for i=1:3
%             if i==1
%                 A{1}=wavelet(A0,kernel{1});
%                 if flag==1
%                    T0=IUWT(A0);
%                    T{1}=IUWT(A{1}); 
%                    W{1}=T0-T{1};
%                 else
%                    W{1}=A0-A{1}; 
%                 end
%             else
%                 A{i}=wavelet(A{i-1},kernel{i});
%                 if flag==1
%                    T{i}=IUWT(A{i}); 
%                    W{i}=T{i-1}-T{i};
%                 else
%                    W{i}=A{i-1}-A{i}; 
%                 end
%             end
%         end 
%         W2(:,:,imgcnt) = W{2};
%         W3(:,:,imgcnt) = W{3};
%     end
    
        A0 = double(img);
        %W1
        A1 = conv3_cr_gpu(A0, kernel{1});

        %W2
        A0 = conv3_cr_gpu(A1, kernel{2});
        W2 = A1-A0;

        %W3
        A1 = conv3_cr_gpu(A0, kernel{3});
        W3 = A0-A1;
        

%     for i=1:3
%         deta=median(median( abs(W{i}-median(median(W{i})).*ones(row,column)) ));
%         t=coef*deta/0.67;
%         W{i}=(W{i}>=t).*W{i};
%     end

%     if nargout==1
%         varargout{1}=W{3};
%     elseif nargout==2
%         varargout{1}=W{2};
%         varargout{2}=W{3};
%     elseif nargout==3
%         varargout{1}=W{2};
%         varargout{2}=W{3};
%         varargout{3}=A{3};
%     else
%         varargout{1}=W{1};
%         varargout{2}=W{2};
%         varargout{3}=W{3};
%         varargout{4}=A{3};
%     end   
end


function A=wavelet(I,kernel)
    % wavelet transform    
    [row column]=size(I);
    extend=(length(kernel)-1)/2;
    I=padarray(I,[extend extend],'symmetric');
    I1=conv2(kernel,kernel',I,'same');
    A=I1(extend+1:extend+row,extend+1:extend+column);
end


function T=IUWT(A)
    % Perform a T-transform
    h=[1/16 1/4 3/8 1/4 1/16];
    t1=sum(h);
    t2=sum(h.^2);
    t3=sum(h.^3);
    c=(7*t2)/(8*t1)-t3/(2*t2);
    b=2*sqrt(t1/t2);
    T=b*sign(A+c).*sqrt(abs(A+c));
end