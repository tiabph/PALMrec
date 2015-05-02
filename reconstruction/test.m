imgsize = [1500,2000];
pointnum = 100000;

tic
timg = zeros(imgsize);
timg = AddGaussian2D(timg, rand(1,pointnum).*imgsize(2), ...
    rand(1,pointnum).*imgsize(1),100+rand(1,pointnum).*100, ...
    rand(1,pointnum)+0.5);
toc
imagesc(timg);
colormap gray