function result = EvalPointSet(trans,p1,p2)
%assum: length(p1)<=length(p2)
    tp1 = GetPointSet(p1,trans);
    len1 = size(p1,1);
    len2 = size(p1,2);
    tlen = round(len1.*0.5);
    dist = zeros(1,len1);
    for m=1:len1
        tdist = (tp1(m,1)-p2(:,1)).^2 + (tp1(m,2)-p2(:,2)).^2;
        dist(m)=min(tdist);
    end
    dist=sort(dist);
    result = sum(dist(1:tlen));
end

function result = GetPointSet(x,trans)
    result = x*[trans(1) trans(2)
                trans(3) trans(4)];
    result(:,1)=result(:,1)+trans(5);
    result(:,2)=result(:,2)+trans(6);
end