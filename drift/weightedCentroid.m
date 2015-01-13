function [x y] = weightedCentroid(I)
    s=size(I);
    if(length(s)<3)
        s(3)=1;
    end
    [xx yy] = meshgrid(1:s(2), 1:s(1));
    x = zeros(s(3),1);
    y = zeros(s(3),1);
    for m=1:s(3)
        timg = double(I(:,:,m));
        x(m) = sum(timg(:).*xx(:))/sum(timg(:));
        y(m) = sum(timg(:).*yy(:))/sum(timg(:));
    end
end