function result = parseArgin(inputbuf, result)
    if(nargin <2)
        result = [];
    end
    
    if(iscell(inputbuf))
        for m=1:length(inputbuf)
            result = parseLine(inputbuf{m}, result);
        end
    elseif(ischar(inputbuf))
            result = parseLine(inputbuf, result);
    end
end

function result = parseLine(inputline, result)
    pos1 = strfind(inputline,'=');
    pos1 = pos1(1);
    key = inputline(1 : pos1-1);
    val = inputline(pos1+1 : end);
    key = strtrim(key);
    val = strtrim(val);
    if(~isnan(str2double(val)))
        val = str2double(val);
    end
    pos = strfind(key,'.');
    if(isempty(pos))
        result = setfield(result,key, val);
    else
        key1 = key(1:pos-1);
        key2 = key(pos+1:end);
        if(isfield(result, key1))
            tstruct = getfield(result,key1);
        else
            tstruct = [];
        end
        tstruct = setfield(tstruct,key2, val);
        result = setfield(result,key1, tstruct);
    end
end