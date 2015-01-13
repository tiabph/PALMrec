function B = CalDrift_GPU(I)
%     [wx wy] = weightedCentroid(I);
    sigma = 2;
    [ret, y x N BG S CRLBy CRLBx CRLBn CRLBb CRLBs LogL] = evalc('mGPUgaussMLE(single(I),sigma,20,2);');
%     x=smooth(x,50,'rlowess');
%     y=smooth(y,50,'rlowess');
%     x=x-x(1);
%     y=y-y(1);
    B = [x(:) y(:)];
end