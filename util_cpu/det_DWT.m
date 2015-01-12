function [W2 W3] = det_DWT(img, flag)
%wavelet tranform for partical detection
%flag for use cup(1) or gpu(0)
    if(nargin<2)
        flag = 1;
    end
    kernel{1}=[1/16,1/4,3/8,1/4,1/16];
    kernel{2}=[1/16,0,1/4,0,3/8,0,1/4,0,1/16];
    kernel{3}=[1/16,0,0,0,1/4,0,0,0,3/8,0,0,0,1/4,0,0,0,1/16];
    
    if(flag)
        [W2 W3] = det_DWT_cpu_pad(double(img),kernel{1}, kernel{2},kernel{3});
    else
        [W2 W3] = det_DWT_gpu(img);
    end
end