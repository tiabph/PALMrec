function WriteINI(data, inifile)
    fp = fopen(inifile,'w');
    names = fieldnames(data);
    for m=1:length(names)
        key = names{m};
        val = getfield(data, key);
        if(~isstruct(val))
            fprintf(fp, '%s=%s\n', key, num2str(val));
        end
    end
    for m=1:length(names)
        key = names{m};
        val = getfield(data, key);
        if(isstruct(val))
            writeSection(fp, key, val)
        end
    end
    fclose(fp);
end

function writeSection(fp, secname, data)
    fprintf(fp, '[%s]\n', secname);
    names = fieldnames(data);
    for m=1:length(names)
        key = names{m};
        val = getfield(data, key);
        fprintf(fp, '%s=%s\n', key, num2str(val));
    end
end