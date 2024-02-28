clear;clc;close all; 
%%%%%%% Set up path %%%%%%%%%%%%%%%%%%
mask_base_dir = 'Z:\Jessica\raw\160623_raw\2B_20x'; %base directory for segmented masks
mask_csv_subdir = '\composites\csv'; %directory name where cellpose mask csv files aare stored under each shot
projectpath='Z:\Jessica\raw\160623_raw\2B_20x'; %base dir of raw microscopy images
imagepath='Z:\Jessica\raw\160623_raw';
experimentpath='\2B_20x\'; %dir under each shot where images are stored
figure_dir = 'Z:\Jessica\tracking_code\160523\figures_1'; %dir for figures generated
data_dir = 'Z:\Jessica\tracking_code\160523\data_1';%directory for matrix data generated
add_path; %add the required dependencies to path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make data and figure directory if not already exist
if ~exist(figure_dir,'dir')
    mkdir(figure_dir);
end

if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

%the shots to be tracked, change this if diff shots are used
shots_row = {'1', '1','3','3'}; shots_col = {'6','6','1', '1'}; shots_sites = {'1','2', '1','2'};
delimiter = {'_'};
%format shots number
shot_cat = strcat(shots_row,delimiter,shots_col,delimiter,shots_sites);


%perform nuclear tracking and match mask from cellpose to tracking results
for index = 1%:numel(shot_cat)
    tic
    shot_num = [str2double(shots_row(index)), str2double(shots_col(index)), str2double(shots_sites(index))];
    file_path = char(fullfile(mask_base_dir, shot_cat(index), mask_csv_subdir));
    
    %nucler tracking
    nuc_tracking = Timelapse20x(shot_num(1), shot_num(2),shot_num(3),projectpath,imagepath,experimentpath, 1,25);
    
    total_frames = length(struct2cell(dir(fullfile(file_path, '*.csv'))));
    %add coordination speed analysis results to trace data
    trace_motility = CoordinationSpeed(nuc_tracking,shot_num(1), shot_num(2),1,total_frames-1);
    save([data_dir,filesep,strcat(shot_cat{index},'_tracking_varied_model')],"trace_motility"); 
    %mask tracking
    mask_track = mask_matching(file_path, trace_motility, shot_cat{index});
    save([data_dir,filesep,strcat(shot_cat{index},'_mask_tracking_varied_model')],"mask_track"); %save track matrix

    %separate out centroid matrix for masks and nuclei
    mask_track_x = mask_track(:,:,size(mask_track,3)-4);
    mask_track_y = mask_track(:,:,size(mask_track,3)-3);
    nuc_track_x = mask_track(:,:,1);
    nuc_track_y = mask_track(:,:,2);
    count_tracks = sum(~isnan(mask_track_x),2);
    count_track_nuc = sum(~isnan(nuc_track_x),2);
    figure;histogram(count_tracks,36)
    xlabel({'Track length', '(frames)'});
    ylabel({'Frequency'});
    figure;histogram(count_track_nuc,36)
    xlabel({'Track length', '(frames)'});
    ylabel({'Frequency'});
    mask_area = mask_track(:,:,size(mask_track,3)-5);
    relative_area = bsxfun(@minus,mask_area,mask_area(:,1));
    cell_velocity = mask_track(:,:,size(mask_track,3)-12);
    %get x,y coordinates relative to the starting point of each track
    relative_x = bsxfun(@minus,mask_track_x,mask_track_x(:,1));
    relative_y = bsxfun(@minus,mask_track_y,mask_track_y(:,1));
    relative_x_nuc = bsxfun(@minus,nuc_track_x,nuc_track_x(:,1));
    relative_y_nuc = bsxfun(@minus,nuc_track_y,nuc_track_y(:,1));
    save([figure_dir,filesep,strcat(shot_cat{index},'_relative_tracks_varied_model')],'relative_x','relative_y','relative_x_nuc', 'relative_y_nuc');
   
    %generate trajectory plot using mask tracks (mask centroid coordinates)
    relative_mask_pos = figure;
    plot(relative_x', relative_y')
    axis tight
    xlabel({'x position', '(pixels)'});
    ylabel({'y position', '(pixels)'});
    saveas(relative_mask_pos,[figure_dir,filesep,strcat(shot_cat{index},'_relative_mask_tracks_varied_model')],'png');

    abs_mask_pos = figure;
    plot(mask_track_x', mask_track_y')
    axis tight
    xlabel({'x position', '(pixels)'});
    ylabel({'y position', '(pixels)'});
    saveas(abs_mask_pos,[figure_dir,filesep,strcat(shot_cat{index},'_abs_mask_tracks_varied_model')],'png');

    %generate trajectory plots (relative and absolute trajectory) based on
    %nucleus centroid coordinates 
    relative_nuc_pos = figure;
    plot(relative_x_nuc', relative_y_nuc')
    axis tight
    xlabel({'x position', '(pixels)'});
    ylabel({'y position', '(pixels)'});
    saveas(relative_nuc_pos,[figure_dir,filesep,strcat(shot_cat{index},'_relative_nuclear_tracks_varied_model')],'png');

    abs_nuc_pos = figure;
    plot(nuc_track_x', nuc_track_y')
    axis tight
    xlabel({'x position', '(pixels)'});
    ylabel({'y position', '(pixels)'});
    saveas(abs_nuc_pos,[figure_dir,filesep,strcat(shot_cat{index}, '_abs_nuclear_tracks_varied_model')],'png');
    toc
end 