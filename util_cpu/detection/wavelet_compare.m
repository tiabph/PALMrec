function A=wavelet_compare(I,kernel)
%% wavelet transform    
[row column]=size(I);
extend=(length(kernel)-1)/2;
I=padarray(I,[extend extend],'symmetric');
I1=conv2(kernel,kernel',I,'same');
A=I1(extend+1:extend+row,extend+1:extend+column);