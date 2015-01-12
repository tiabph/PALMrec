mex det_DWT_cpu.cpp
t=whos('img');
if(length(t)==0)
    img = tiffread('e:\a647.tif');
end

kernel{1}=[1/16,1/4,3/8,1/4,1/16];
kernel{2}=[1/16,0,1/4,0,3/8,0,1/4,0,1/16];
kernel{3}=[1/16,0,0,0,1/4,0,0,0,3/8,0,0,0,1/4,0,0,0,1/16];
wt1 = kernel{1}'*kernel{1};
wt2 = kernel{2}'*kernel{2};
wt3 = kernel{3}'*kernel{3};

tic
[W2 W3] = det_DWT_cpu(double(img), wt1, wt2, wt3);
toc

tic
[tW2 tW3] = det_DWT(img, 0);
toc