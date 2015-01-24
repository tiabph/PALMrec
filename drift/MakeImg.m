function img = MakeImg(pos, simg, zoom)
    pos = pos .*zoom;
    simg = simg * zoom;
    
    img = zeros(simg);
    for m=1:size(pos,1)
        tx = round(pos(m,1));
        ty = round(pos(m,2));
        if(tx >0 && tx <=simg(2) && ty>0 && ty <=simg(1))
            img(ty,tx) = img(ty,tx) +1;
        end
    end
end