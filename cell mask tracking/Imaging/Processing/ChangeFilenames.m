%% Change filenames from IXµ raw format to a file format to be used with Min's
% tracking code


for rows=2:7
    for cols=2:6
        for sites=1
            FormatFiles_1site_3color(rows,cols,sites);
            %FormatFiles(rows,cols,sites);
        end
    end
end

disp('done!');