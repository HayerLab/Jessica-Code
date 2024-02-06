%% Change filenames from Axon raw format (acquired using Matlab)
% to a file format to be used with Min's tracking code, e.g.
% 2_1_1_DAPI_1.tif (row_col_site_DAPI_frame.tif).
% 10 July 2014 Arnold
clear; clc;
rawdir='D:\data\140708_Axon\Mitofission\';
targetdir='D:\data\140708_Axon\Mitofission\';

for rows=3:7
    alp=char(64+rows);
    for cols=1:12
        for sites=1
            if cols<10
            tempwell=[alp,'0',num2str(cols)];% alp= alphabetic letter, '0' , num2str fxn assigns a character to the value
            else tempwell=[alp,num2str(cols)];
            end
            DAPIfiles=dir([rawdir,'\',tempwell,'\','*DAPI*']);
            for frames=1:length(DAPIfiles)
                DAPIOld=[rawdir,tempwell,'\',DAPIfiles(frames).name];
                DAPINew=[targetdir,num2str(rows),'_',num2str(cols),'_',num2str(sites),'_DAPI_',num2str(frames),'.tif'];
                system(['move "',DAPIOld,'" "',DAPINew,'"']);
                disp([num2str(rows),' ',num2str(cols),' ',num2str(frames)]);
                %movefile(DAPIOld,DAPINew,'f');
            end
        end
    end
end
disp('done!');