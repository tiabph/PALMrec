function param = CheckParam(param)
% check and fill empty parts in param
% need full path in param 
    if CheckFieldEmpty(param, 'fullpath')
        return
    else
        fullpath = param.fullpath;
        [pathstr, name, ext] = fileparts(fullpath);
        %file info
        param.filename = [name ext];
        param.filepath = [pathstr '\'];
        param.name = name;
        param.ext = ext;
        
        %output file info
        if CheckFieldEmpty(param, 'outpath')
            param.outpath = param.filepath;
        end
        if CheckFieldEmpty(param, 'outfile')
            param.outfile = [name param.outfile_postfix ext];
        end
        
        %drift
        if CheckFieldEmpty(param.drift, 'filepath')
            param.drift.filepath = [param.filepath name param.drift.postfix ext];
        end
    end
end

function ret = CheckFieldEmpty(param, filedname)
    if ~isfield(param, filedname) || isempty(getfield(param, filedname))
        ret = 1;
    else
        ret = 0;
    end
end