rows=1
cols=[3 4 5]  %[3 7]
sites=1%:2%:8 %1:9
% %%
% shot={}; 
% i=0;
% for rows=2:3
%     for cols=2:11
%         i=i+1
%         shot{i}=[num2str(rows),' ',num2str(cols),' 1;'];
%         shot=shot';
%     end 
% end
%%
manualwells = [
    1 1 1;
    1 2 1;
    1 3 1;
    2 5 1;
    2 6 1;
    2 7 1;
    2 8 1;
    2 9 1;
    2 10 1;
    3 9 1;
    3 10 1;
    3 11 1;
    3 12 1;
    4 1 1;
    4 2 1;
    5 1 1;
    6 1 1;
    6 2 1;
    6 3 1;
    6 4 1;
    6 5 1;
    6 6 1;
    6 7 1;
    6 8 1;
    6 9 1;
    6 10 1;
    6 11 1;
    6 12 1;
    7 1 1;
    7 2 1;
    7 3 1;
    7 4 1;
    7 5 1;
    7 6 1;
    7 7 1;
    7 8 1;
    7 9 1;
    7 10 1;
    7 11 1;
    7 12 1;
    8 1 1;
    8 2 1;
    8 3 1;
    8 4 1;
    8 5 1;
    8 6 1;
    8 7 1;
    8 8 1;
    8 9 1;
    8 10 1;
    8 11 1;
    8 12 1;
    ];



manualcontrol=0;
numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
shots=numrows*numcols*numsites;


if manualcontrol==1
    shots=size(manualwells,1);
end
time1=tic;
for shot=1:shots
    if manualcontrol==1
        row=manualwells(shot,1);
        col=manualwells(shot,2);
        site=manualwells(shot,3);
    else
        siteidx=mod(shot,numsites);
        if siteidx==0
            siteidx=numsites;
        end
        site=sites(siteidx);
        colidx=mod(ceil(shot/numsites),numcols);
        if colidx==0
            colidx=numcols;
        end
        col=cols(colidx);
        rowidx=ceil(shot/(numcols*numsites));
        row=rows(rowidx);
    end
    fprintf([num2str(row),'_',num2str(col),'_',num2str(site),'\n']);
    
    %%% Format Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FormatFiles(row,col,site);
    
    %%% Timelapse %%%%%%%%%%%%%%%%%%%%%%%%%
    Timelapse4x(row,col,site);
    
    %%% Immunostain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Immunostain_1_CalcBleedthroughRate(row,col,site);
    %Immunostain(row,col,site);
    %Immunostain_2_AddToTimelapse(row,col,site);
end
toc(time1)