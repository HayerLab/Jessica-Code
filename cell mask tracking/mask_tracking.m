clear;clc;close all; 
%%%%%%% Set up path %%%%%%%%%%%%%%%%%%
mask_base_dir = 'Z:\Jessica\segment'; %base directory for segmented masks
mask_csv_subdir = '\Composites\aligned\tif\csv'; %directory name where cellpose mask csv files aare stored under each shot
projectpath='Z:\Jessica\raw'; %base dir of raw microscopy images
imagepath='Z:\Jessica\raw';
experimentpath='\DAPI_YFP\'; %dir under each shot where images are stored
figure_dir = 'Z:\Jessica\tracking_code\figures'; %dir for figures generated
data_dir = 'Z:\Jessica\tracking_code\data';%directory for matrix data generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make data and figure directory if not already exist
if ~exist(figure_dir,'dir')
    mkdir(figure_dir);
end

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

%the shots to be tracked, change this if diff shots are used
shots_row = [3]; shots_col = [8]; shots_sites = [1 2]; 
shot = {};
shot_cat = {};
ind = 1;

%format shots number so can directly input into tracking function
for i = 1:numel(shots_row)
    for j = 1:numel(shots_col)
        for k = 1:numel(shots_sites)
            shot{ind} = [shots_row(i), shots_col(j), shots_sites(k)];
            shot_str = strcat(num2str(shots_row(i)),'_',num2str(shots_col(j)),'_',num2str(shots_sites(k)));
            shot_cat{ind} =shot_str;
            ind = ind + 1;
        end
    end
end


%perform nuclear tracking and match mask from cellpose to tracking results
for index = 1:numel(shot)
    tic
    shot_num = shot{index};
    file_path = char(fullfile(mask_base_dir, shot_cat(index), mask_csv_subdir));
    
    %nucler tracking
    nuc_tracking = Timelapse20x(shot_num(1), shot_num(2),shot_num(3),projectpath,imagepath,experimentpath);

    total_frames = length(struct2cell(dir(fullfile(file_path, '*.csv'))));
    %add coordination speed analysis results to trace data
    trace_motility = CoordinationSpeed(nuc_tracking,shot_num(1), shot_num(2),1,total_frames-1);
    
    %mask tracking
    mask_track = mask_matching(file_path, trace_motility, shot_cat{index});
    save(char(strcat(data_dir,shot_cat{index},'_mask_tracking')),"mask_track"); %save track matrix

    %separate out centroid matrix for masks and nuclei
    mask_track_x = mask_track(:,:,size(mask_track,3)-5);
    mask_track_y = mask_track(:,:,size(mask_track,3)-4);
    nuc_track_x = mask_track(:,:,1);
    nuc_track_y = mask_track(:,:,2);
    
    %get x,y coordinates relative to the starting point of each track
    relative_x = bsxfun(@minus,mask_track_x,mask_track_x(:,1));
    relative_y = bsxfun(@minus,mask_track_y,mask_track_y(:,1));
    relative_x_nuc = bsxfun(@minus,nuc_track_x,nuc_track_x(:,1));
    relative_y_nuc = bsxfun(@minus,nuc_track_y,nuc_track_y(:,1));
    save(char(figure_dir,strcat(shot_cat{index},'_relative_tracks')),'relative_x','relative_y','relative_x_nuc', 'relative_y_nuc');
   
    %generate trajectory plot using mask tracks (mask centroid coordinates)
    relative_mask_pos = figure;
    plot(relative_x', relative_y')
    xlabel({'x position', '(pixels)'});
    ylabel({'y position', '(pixels)'});
    saveas(relative_mask_pos,char(figure_dir,strcat(shot_cat{index},'_relative_mask_tracks')),'png');

    abs_mask_pos = figure;
    plot(mask_track_x', mask_track_y')
    xlabel({'x position', '(pixels)'});
    ylabel({'y position', '(pixels)'});
    saveas(abs_mask_pos,char(strcat(figure_dir,shot_cat{index},'_abs_mask_tracks')),'png');

    %generate trajectory plots (relative and absolute trajectory) based on
    %nucleus centroid coordinates 
    relative_nuc_pos = figure;
    plot(relative_x_nuc', relative_y_nuc')
    xlabel({'x position', '(pixels)'});
    ylabel({'y position', '(pixels)'});
    saveas(relative_nuc_pos,char(strcat(figure_dir,shot_cat{index},'_relative_nuclear_tracks')),'png');

    abs_nuc_pos = figure;
    plot(nuc_track_x', nuc_track_y')
    xlabel({'x position', '(pixels)'});
    ylabel({'y position', '(pixels)'});
    saveas(abs_nuc_pos,char(strcat(figure_dir,shot_cat{index}, '_abs_nuclear_tracks')),'png');
    toc
end 