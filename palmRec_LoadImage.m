function databuf = palmRec_LoadImage(databuf, param)
    if(~isfield(param, 'filepath') || ~isfield(param, 'filename'))
        disp('palmRec_LoadImage need param has fields below:');
        disp('param.filepath');
        disp('param.filename');
    end
    databuf.path = param.filepath;
    databuf.imgfile = param.filename;
    databuf.img = single(tiffread([databuf.path databuf.imgfile]));
    databuf.height = size(databuf.img, 1);
    databuf.width = size(databuf.img, 2);
    databuf.imglen = size(databuf.img, 3);
end