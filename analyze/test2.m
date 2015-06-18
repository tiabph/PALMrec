%% localization precision analysis test ver2

idxlist = find(databuf.linkData_sp(:,2) == databuf.imglen);

pointnum = length(idxlist);
pbuf = cell(pointnum, 1);
for m=1:pointnum
    ret = GetLinkChain(databuf.linkData_chain, idxlist(m));
    pbuf{m} = databuf.fitInfo(ret,[1 2 7]);
end

frameBuf = cell(databuf.imglen, 1);
for m=1:databuf.imglen
    temp = zeros(pointnum, 3);
    for n=1:pointnum
        temp(n,:) = pbuf{n}(m,:);
    end
    frameBuf{m} = temp;
end

%% registration
basePoints = frameBuf{1}(:,1:2);
regBuf = cell(databuf.imglen, 1);
mask = true(pointnum, 1);

for m=1:databuf.imglen
    temp = frameBuf{m}(:, 1:2);
    [tresult ttransmat tmask] = RegPointsPair(basePoints, temp);
    
    regBuf{m} = tresult;
    mask = mask & tmask;
end

regPointBuf = cell(pointnum, 1);
for m=1:pointnum
    temp = zeros(databuf.imglen, 3);
    for n=1:databuf.imglen
        temp(n,1:2) = regBuf{n}(m,:);
    end
    temp(:, 3) = pbuf{m}(:,3);
    regPointBuf{m} = temp;
end

%% calculate localization precision
stdbuf = zeros(pointnum, 2);
for m=1:pointnum
    temp = zeros(databuf.imglen, 3);
    for n=1:databuf.imglen
        temp(n,1:2) = regBuf{n}(m,:);
    end
    temp(:, 3) = pbuf{m}(:,3);
    
    [muhatx,sigmahatx] = normfit(temp(:,1));
    [muhaty,sigmahaty] = normfit(temp(:,2));
    stdbuf(m,1) = sigmahatx;
    stdbuf(m,2) = sigmahaty;
end

std_max = max(stdbuf,[],2);
[std_sorted std_idxlist] = sort(std_max);

%% plot result
figure(1)
imagesc(databuf.img(:,:,1));
colormap gray
hold on
for m=1:databuf.imglen
    temp = regBuf{m}(mask, 1:2);
    plot(temp(:,1), temp(:,2),'r.');
end

hold off

figure(2)
for m=1:9
    subplot(3,3,m);
    plot((regPointBuf{std_idxlist(m)}(:,1) - mean(regPointBuf{std_idxlist(m)}(:,1)))*param.fitting.pixelsize, ...
        (regPointBuf{std_idxlist(m)}(:,2)-mean(regPointBuf{std_idxlist(m)}(:,2)))*param.fitting.pixelsize, '.');
    title(['std: x' num2str(stdbuf(std_idxlist(m),1)*param.fitting.pixelsize) ' y' num2str(stdbuf(std_idxlist(m),2)*param.fitting.pixelsize)]);
end