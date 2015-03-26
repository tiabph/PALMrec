function result = ReadINI(inifile, autoconverse, result)
    if(nargin<2)
        autoconverse = 1;
    end
    if(nargin<3)
        result = [];
    end
    fp = fopen(inifile,'r');
    while(1)
        templine = fgetl(fp);
        if(length(templine)==1 && templine==-1)
            break;
        end
        [rtype key value] = parseLine(templine, autoconverse);
        if(rtype ==1)
            result = setfield(result, key, value);
        end
        if(rtype ==2)
            sectionname = key;
            tstruct = [];
            while(1)
                templine = fgetl(fp);
                if(length(templine)==1 && templine==-1)
                    if(~isempty(tstruct))
                        result = setfield(result, sectionname, tstruct);
                    end
                    break;
                end
                [rtype key value] = parseLine(templine, autoconverse);
                if(rtype ==1)
                    tstruct = setfield(tstruct, key, value);
                end
                if(rtype ==2)
                    result = setfield(result, sectionname, tstruct);
                    sectionname = key;
                    tstruct = [];
                end
            end
        end
    end
    
    fclose(fp);
end

%type: 0:comment, 1:kay-value pair, 2:section
function [rtype key value] = parseLine(line, convflag)
    line = strtrim(line);
    if(line(1) == ';' || line(1) == '#')
        rtype = 0;
        key = [];
        value = [];
        return
    elseif(line(1) == '[' && line(end) == ']')
        rtype = 2;
        key = line(2:end-1);
        value = [];
        return
    end
    rtype = 1;
    pos = strfind(line, '=');
    pos = pos(1);
    str_key = line(1:pos-1);
    str_val = line(pos+1:end);
    key = strtrim(str_key);
    value = strtrim(str_val);
    if(convflag>0 && ~isnan(str2double(value)))
        value = str2double(value);
    end
end