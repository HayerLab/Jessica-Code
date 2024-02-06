rows=2:7
cols=1:12    %[3 7]
sites=1 %1:9

manualwells = [
    2 2 1;
    2 3 1;
    2 4 1;
    2 5 1;
    2 6 1;
    2 7 1;
    2 9 1;
    2 10 1;
    2 11 1;
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
parfor shot=1:shots
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
    Timelapse_axon(row,col,site);
    
    %%% Immunostain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Immunostain_1_CalcBleedthroughRate(row,col,site);
    %Immunostain(row,col,site);
    %Immunostain_2_AddToTimelapse(row,col,site);
end
toc(time1)