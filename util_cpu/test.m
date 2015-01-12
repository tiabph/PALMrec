%test scrip
t=whos('img');
if(length(t)==0)
    img = tiffread('e:\a647.tif');
end
% img = img(:,:,1:1000);

%%
[detectResults W2 tW2] = DWTParticalDetection(img, []);

%%
points = zeros(0,3);
for m=1:length(detectResults)
    points = cat(1, points, detectResults{m});
end

figure(1)
plot(points(:,1), points(:,2), '.');
imgsize = size(img);
zoom = 8;
imgbuf = zeros(imgsize(1)*zoom, imgsize(2)*zoom);
for m=1:size(points,1)
    imgbuf(floor(points(m,2)*zoom), floor(points(m,1)*zoom)) = imgbuf(floor(points(m,2)*zoom), floor(points(m,1)*zoom)) +1;
end
figure(2)
imagesc(imgbuf)
colormap gray