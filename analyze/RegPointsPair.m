function [result transmat mask] = RegPointsPair(points_base, points_test)
    sigma_threshold = 4;
    regtype = 0;
    
    pointnum = size(points_base, 1);
    distance = points_base - points_test;
    [muhatx,sigmahatx] = normfit(distance(:,1));
    [muhaty,sigmahaty] = normfit(distance(:,2));
%     mask = ones(pointnum, 1);
    maskx = abs(distance(:,1) - muhatx) <= (sigmahatx*sigma_threshold);
    masky = abs(distance(:,2) - muhaty) <= (sigmahaty*sigma_threshold);
    mask = maskx & masky;
    
    pbase = [points_base(mask, 1:2) ones(sum(mask),1)];
    ptest = [points_test(mask, 1:2) ones(sum(mask),1)];
    
    if regtype == 1
        %calculate tranform matrix
        transmat = (ptest' * ptest) \ (ptest' * pbase);

        result = [points_test(:,1:2) ones(pointnum,1)] * transmat;
        result = result(:,1:2);
    else
        result = points_test(:,1:2);
        result(:,1) = result(:,1) + muhatx;
        result(:,2) = result(:,2) + muhaty;
        
        transmat = [1 0 0; 0 1 0; muhatx muhaty 1];
    end
end