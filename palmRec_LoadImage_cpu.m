function databuf = palmRec_LoadImage_cpu(databuf, param)
    if(~isfield(param, 'filepath') || ~isfield(param, 'filename'))
        disp('palmRec_LoadImage need param has fields below:');
        disp('param.filepath');
        disp('param.filename');
    end
    databuf.path = param.filepath;
    databuf.imgfile = param.filename;
    databuf.img = single(TIFFloadFrame_16bit_CPU([databuf.path databuf.imgfile],[1 -1]));
    databuf.height = size(databuf.img, 1);
    databuf.width = size(databuf.img, 2);
    databuf.imglen = size(databuf.img, 3);
end