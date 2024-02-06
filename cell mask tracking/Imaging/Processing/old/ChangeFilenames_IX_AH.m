%% Change filenames from IXµ raw format to a file format to be used with Min's
% tracking code

clear; clc;
for row=1:8
    for col=1:12
        for site=1:3
            %FormatFiles_old1(row,col,site);
            FormatFiles_old(row,col,site);
        end
    end
end