clc
clear
close all
%%%%%%% 
% This script relabels each cell mask with its corresponding cell number.
% Relabeled matrix for each frame is saved. 

mask_base_dir = 'Z:\Jessica\segment\';
mask_csv_subdir = '\Composites\aligned\tif\csv';
data_dir = 'Z:\Jessica\tracking_code\data\varied_model';%directory for matrix data generated
label_dir = 'Z:\Jessica\tracking_code\label_matrix';
shots_row = {'3'}; shots_col = {'8'}; shots_sites = {'1','2'};
delimiter = {'_'};
%format shots number
shot_cat = strcat(shots_row,delimiter,shots_col,delimiter,shots_sites);

%for index = 2
for index = 1:numel(shot_cat)
    mask_track = importdata([data_dir, filesep,strcat(shot_cat{index},'_mask_tracking_varied_model.mat')]);
    mask_label = mask_track(:,:,13);
    disp(shot_cat{index})

    %read cellpose csv files
    file_path = char(fullfile(mask_base_dir, shot_cat(index), mask_csv_subdir));
    files = dir(fullfile(file_path, '*.csv'));
    files = natsortfiles(files); %sort files in natural order
    files = struct2cell(files);
    files = files(1,:);
    %label_mat_first = readmatrix(strcat(file_path, filesep, files{1}));

    for idx = 1:length(files)
        %read mask csv file to retrieve label matrix for cell masks
        mask = readmatrix(strcat(file_path, filesep, files{idx}));
        cell_label_mat = zeros([size(mask,1), size(mask,2)]);

        %relabel cell matrix with cell numbers
        cell_label_mat = cell_number(mask, idx, cell_label_mat, mask_track,1);
        disp(strcat('Frame number:',num2str(idx)))
        save([label_dir,filesep,shot_cat{index},filesep,strcat(shot_cat{index},'_',num2str(idx),'_cell_num_matrix')],"cell_label_mat");
    end
end