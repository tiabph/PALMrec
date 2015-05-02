function img = AddGaussian2D(img, x, y, height, std)
    s = size(img);
    patchsize = 3;
    [xx yy] = meshgrid(1:s(2), 1:s(1));
    
    for m=1:length(x)
        x0=x(m);
        y0=y(m);
        psize = floor(patchsize*std(m));
        xs = max(1, round(x0) - psize);
        xe = min(s(2), round(x0) + psize);
        ys = max(1, round(y0) - psize);
        ye = min(s(1), round(y0) + psize);
        
        result = fun_Gaussian(xx(ys:ye, xs:xe),yy(ys:ye, xs:xe),x0,y0,std(m));
        img(ys:ye, xs:xe) = img(ys:ye, xs:xe) + result.*height(m);
    end
end

function result = fun_Gaussian(x,y,x0,y0,std)
    result = exp(-((x-x0).^2 + (y-y0).^2)/(2*std^2))/(2*pi*std^2);
end