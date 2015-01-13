function B = CalDrift_weightedcentroid(I)
    [x y] = weightedCentroid(I);
%     x=smooth(x,50,'rlowess');
%     y=smooth(y,50,'rlowess');
%     x=x-x(1);
%     y=y-y(1);
    B = [x(:) y(:)];
end