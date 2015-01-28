load drift_dataset.mat

buglen = 10000;
[xx yy] = meshgrid(-20:20, -20:20);
dnum = floor(size(pos,1)/buglen);
dresult = zeros(dnum,2);
imgbuf = cell(dnum);
for m=1:dnum-1
    m
    img1 = MakeImg(pos(((1-1)*buglen+1):(1*buglen),:), simg, 32);
    img2 = MakeImg(pos((m*buglen +1):((m+1)*buglen),:), simg, 32);
    reimg = calCrossCorr(img1, img2, xx,yy);
    B = CalDrift_weightedcentroid(reimg);
    B = B - 10;
    dresult(m+1,:) = B;
    imgbuf{m} = reimg;
end


dresult = dresult - 11;
result_i = zeros(size(dresult));
for m=2:size(dresult,1)
result_i(m,1) = dresult(m,1) + result_i(m-1,1);
result_i(m,2) = dresult(m,2) + result_i(m-1,2);
end

figure(2)
plot(dresult)
