function Wx = det_Thresh(Wx, threshold)
    row = size(Wx, 1);
    column = size(Wx, 2);
    
    t=reshape(Wx,row*column,size(Wx,3));
    mlist = mean(t,1);
    stdlist = std(t,1);
    for i=1:size(Wx,3)
        timg=Wx(:,:,i);
        %replace median with mean, factor of log(2)
        %deta_2=mean( abs(timg(:)-mean(timg(:))));
%         deta=median(abs( timg(:)-median(timg(:)) ));
%         deta = detalist(i);
%         t=threshold*deta/0.67;
        tsh = mlist(i) + stdlist(i)*threshold;
        timg(timg<tsh)=tsh;
        Wx(:,:,i) = timg;
    end
end