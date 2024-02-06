% plate_persistence_screen.m
%     The input datastructure tracedata.mat, but extended by additional layers in the 3rd dimension for
%     motility parameters, for each cell, and each timepoint: 
%     tracedata(:,:,9)  - angle of movement
%     tracedata(:,:,10) - single cell velocity
%     tracedata(:,:,11) - coordination with neighboring cells in the front, located within a 60° sector ahead
%     tracedata(:,:,12) - coordination with lateral neighbors, located within lateral 120° sectors)
%     tracedata(:,:,13) - coordination with neighboring cells in the back, located within a 60° sector behind
%     tracedata(:,:,14) - coordination with all neighboring cells
%     tracedata(:,:,15) - persistence. Col1 : Net distance. Col2 : Total distance. Col3: Net/total distance  
%
% Calculate the ratio of net displacement over total path length. 
% Total path length is the sum of single cell velocity per time point over
% the total Euklidian distance between the first and the last point. 
%
%% Set up
clear;clc; warning off;
SF=1;
EF=25;
root='K:\160831_IXL';
%platenames={'1A','1B','2A','2B','1C','1D','2C','2D','1E','1F','2E','2F'};


%%
%for plate=1:length(platenames)
    sourcedir=[root,filesep,'Results',filesep,'coordination_analysis_160901'];
    targetdir=[root,filesep,'Results',filesep,'CoordinationPersistence'];
    if ~exist(targetdir)
        mkdir(targetdir);
    end
    
    for rows=1:8
        for cols=1:12
            shot=[num2str(rows),'_',num2str(cols),'_1'];
            load([sourcedir,filesep,shot,'_coordination_analysis100µm_1-24.mat']);
            disp(shot);
            % Net distance traveled between SF and EF
            LocStart=[tracedata(:,SF,1) tracedata(:,SF,2)]; % x-y coordinates at time point SF
            LocEnd=[tracedata(:,EF,1) tracedata(:,EF,2)]; % x-y coordinates at time point EF
            StartEnd=LocEnd-LocStart;
            [TotAng TotDist]=cart2pol(StartEnd(:,1),StartEnd(:,2));
            tracedata(:,1,15)=TotDist;

            % Sum of displacement per frame:
            velocity=tracedata(:,:,10);
            tracedata(:,2,15)=sum(velocity(:,1:(EF-1)),2);
            
            % Persistence = Net distance / total distance
            tracedata(:,3,15)=tracedata(:,1,15)./tracedata(:,2,15);
            save([targetdir,filesep,shot,'_coor_100µm_1-24_pers_1-25.mat'],'tracedata');
        end
    end
%end % plate
disp('done!');

