clc
clear
close all
%%%%%%% 
% This script uses the number of occurences of numeric cell label as an
% approximation for cell tracking accuracy.
mask_base_dir = 'Z:\Jessica\segment\';
mask_csv_subdir = '\Composites\aligned\tif\csv';
shots_row = 3; shots_col = 8; shots_sites = [1 2];
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

load('data\mask_track_3-8-1.mat');
mask_label = mask_track(:,:,5);
%count how many columns each cell label occurs in
[unique_num, unique_count] = count_unique(mask_label); 
unique_num = unique_num(2:length(unique_num)); %remove 0, 0 is non-cell mask 
unique_count = unique_count(2:length(unique_count)); % remove # occurence for 0
idx = find(unique_count > 90);
percentage = {};
cell_shift_threshold = 3;

for index = 1
%for index = 1:numel(shot)
    shot_num = shot{index};
    disp(num2str(shot_num))
    file_path = char(fullfile(mask_base_dir, shot_cat(index), mask_csv_subdir));
    files = dir(fullfile(file_path, '*.csv'));
    files = struct2cell(files);
    files = files(1,:);
    label_mat_first = readmatrix(strcat(file_path, filesep, files{i}));
    cell_num_label = zeros(size(label_mat_first,1), size(label_mat_first,2),180);
    %cell_num_label = cell(size(label_mat_first,1), size(label_mat_first,2));
    %cell_num_label = struct([]);
    %tracking_matrix = importdata(char(strcat(file_path,filesep,'cellpose_mask_tracking_',char(shot_cat{index}),'.mat')));
    tracking_matrix = importdata('data\mask_track_3-8-1.mat');
    for i = 1:length(files)
        %read mask csv file to retrieve label matrix for cell masks
        mask = readmatrix(strcat(file_path, filesep, files{i}));
        mask(1, :) = [];
        mask(:,1) = [];
        cell_num_label = cell_number(mask, i, 180, cell_num_label, tracking_matrix);
        %timeit(cell_num_label)
        disp(strcat('Frame number:',num2str(i)))
    end
    %find num of pixels whose corresponding # of cell num does not exceed a
    %threshold value
    % transforming so each row is a slice in row form
    unique_count = zeros(size(cell_num_label,1),size(cell_num_label,2));
    %reshaped_mat = reshape(cell_num_label, [], size(cell_num_label, 3))';
    %unique_count = accumarray(repmat(1:size(reshaped_mat,1),1,size(reshaped_mat,2)).',reshaped_mat(:),[],@(x) numel(unique(x)));
    for i=1:size(cell_num_label,1)
        for j = 1:size(cell_num_label,2)
            unique_count(i,j) = size(unique(cell_num_label(i,j,:)),1);
        end
    end

    threshold_pix = find(unique_count<30 & unique_count>0);
    area = find(unique_count > 1);
    %threshold_pix = all(cellfun(@(cell) le(cell_shift_threshold, unique(cell))&& gt(0, unique(cell)) , cell_num_label(1:end)));
    percentage{end+1} = size(threshold_pix,1)/size(area,1);
end 



