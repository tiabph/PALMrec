function data_w=Detection(A, threshold, coefficient)
A=single(A);

% threshold=parameter.detection.wavelet_threshold;
% coefficient=parameter.detection.wavelet_coefficient;
[W2 W3]=waveletTransform(A,1,threshold);
if coefficient==2
    PJ=W2+W3;
elseif coefficient==1
    PJ=sqrt(W2.*W3);
elseif coefficient==3
    PJ=W2;
else
    PJ=W3;   
end
h=[1/16 1/4 3/8 1/4 1/16];
t1=sum(h);
t2=sum(h.^2);
t3=sum(h.^3);
c=(7*t2)/(8*t1)-t3/(2*t2);
b=2*sqrt(t1/t2);   
PJ1=(PJ.^2)+(b^(-2))-c;
PJ1=(PJ1>0).*PJ1; 
data_w=PJ1;  
