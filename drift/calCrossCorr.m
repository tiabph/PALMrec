function reimg = calCrossCorr(img1, img2, xx,yy)
    reimg = zeros(size(xx));
    s1 = size(img1);
    s2 = size(img2);
    for m=1:length(xx(:))
        dx = xx(m);
        dy = yy(m);
        
        sx1 = max(1+dx,1);
        ex1 = min(s1(2)+dx, s1(2)); 
        sy1 = max(1+dy,1);
        ey1 = min(s1(1)+dy, s1(1)); 
        sx2 = max(1-dx,1);
        ex2 = min(s2(2)-dx, s2(2)); 
        sy2 = max(1-dy,1);
        ey2 = min(s2(1)-dy, s2(1)); 
        
        timg = img1(sy1:ey1,sx1:ex1) .* img2(sy2:ey2,sx2:ex2);
        reimg(m) = mean(timg(:));
    end
end