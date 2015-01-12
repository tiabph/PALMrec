function databuf = palmRec_Save(databuf, param)
    path = param.filepath;
    filename = param.filename;
    amp = param.reconstruction.amp;
    filename = [filename(1:end-4) '_' param.reconstruction.recType '_' num2str(amp) 'X.tif'];
    databuf.outpath = path;
    databuf.outfile = filename;
    tiffwrite(databuf.recimg, [path filename]);
end