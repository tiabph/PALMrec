function [order] = FindPointOrder(p1,p2,MaxDist)
    l1 = size(p1,1);
    l2 = size(p2,1);
    resultlen = min(l1,l2);
    order = zeros(resultlen,2);
    orderlen = 0;
    DistMatrix = zeros(l1,l2);
    for n=1:l1
        for m=1:l2
            DistMatrix(n,m) = sum((p1(n,:)-p2(m,:)).^2);
        end
    end
    for n=1:resultlen
        [result,remmax,DistMatrix] = FindPointSet(DistMatrix);
        orderlen = orderlen+1;
        order(orderlen,:) = result;
        if(remmax >= MaxDist)
            break;
        end
    end
    order=order(1:orderlen,:);
end

function [result,remmax,newmatrix] = FindPointSet(DistMatrix)
    mindist = min(min(DistMatrix));
    [r,c]=find(DistMatrix == mindist);
    r=r(1);
    c=c(1);
    result = [r c];
    newmatrix = DistMatrix;
    newmatrix(r,:) = nan;
    newmatrix(:,c) = nan;
    remmax = min(min(newmatrix));
end