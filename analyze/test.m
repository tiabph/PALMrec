%% localization precision analysis

pointbuf = databuf.fitInfo;
pointbuf = pointbuf(:, [1 2 7 9]);
framenum = databuf.imglen;
pbuf = cell(framenum,1);
spbuf = zeros(framenum,1);
spbuf(1)=1;
for m=2:framenum
    spbuf(m) = find(pointbuf(:,4)==m, 1);
    pbuf{m-1} = pointbuf(spbuf(m-1):(spbuf(m)-1),:); 
    m
end
pbuf{end} = pointbuf(spbuf(end):end, :); 

%% registration

base_pointset = pbuf{1}(:,1:2);
pbuf_reg = cell(framenum,1);
parfor m=1:framenum
    tempbuf = pbuf{m}(:,1:2);
    [tx1 result flag order tfval] = GetPointsReg(tempbuf,base_pointset);
    pbuf_reg{m} = tx1;
    m
end

%% display result
figure(1)
hold off
for m=1:framenum
    temp = pbuf{m}(:,1:2);
    plot(temp(:,1),temp(:,2),'.');
    hold on
end
hold off

figure(2)
hold off
for m=1:framenum
    temp = pbuf_reg{m}(:,1:2);
    plot(temp(:,1),temp(:,2),'.');
    hold on
end
hold off