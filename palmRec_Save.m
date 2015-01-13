function databuf = palmRec_Save(databuf, param)
    path = param.filepath;
    filename = param.filename;
    amp = param.reconstruction.amp;
    
    str_dc = '';
    if isfield('DriftCorrectionFlag', databuf) && DriftCorrectionFlag>0
        str_dc = '_DC';
    end
    filename = [filename(1:end-4) '_' param.reconstruction.recType '_' num2str(amp) 'X' str_dc '.tif'];
    databuf.outpath = path;
    databuf.outfile = filename;
    tiffwrite(databuf.recimg, [path filename]);
end