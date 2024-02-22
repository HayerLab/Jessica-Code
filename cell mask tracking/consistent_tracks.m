clc; clear; close all;
% this files removes any cells that are not presetn in all frames (i.e.
% only cells with track length of 180 are preserved)

%%%%%% set up paths %%%%%%%
label_dir = 'Z:\Jessica\tracking_code\label_matrix\single_nuclei';
label_im_dir = 'Z:\Jessica\tracking_code\label_matrix\single_nuclei\images';
base_dir = 'Z:\Jessica\segment'; %base directory for segmented masks
image_subdir = '\Composites\aligned\tif\';
rawdir = 'Z:\Jessica\raw\160519_CDH5_mCitrine_Hoechst_1min_interval_20x';
shots_row = {'3'}; shots_col = {'8'}; shots_sites = {'1','2'};
delimiter = {'_'};
%format shots number
shot_cat = strcat(shots_row,delimiter,shots_col,delimiter,shots_sites);

for index = 1:numel(shot_cat)
    %path where cell matrix is stored
    matrix_storage_folder = [label_dir,filesep,'data',filesep,shot_cat{index},filesep,'matrix'];
    
    consistent_label_folder = [label_dir,filesep,'data',filesep,shot_cat{index},filesep,'consistent'];
    if ~exist(consistent_label_folder,'dir')
        mkdir(consistent_label_folder);
    end

    consistent_im_folder = [consistent_label_folder,filesep,'images'];
    if ~exist(consistent_im_folder,'dir')
        mkdir(consistent_im_folder);
    end

    %imported filtered labels list
    filtered_labels_list = importdata([label_dir,filesep,'data',filesep, shot_cat{index},filesep,'labels',filesep,'filtered_labels_no_clusters_',shot_cat{index},'.mat']);
    tic

    %find all cell labels that is present across all frames
    all_labels = struct2cell(filtered_labels_list);
    max_size = max(cellfun(@numel,all_labels));
    padded = cellfun(@(x) [x; zeros(max_size - numel(x), 1)], all_labels, 'un', 0);
    all_labels = cell2mat(padded);
    unique_val = unique(all_labels);
    all_labels = squeeze(all_labels)';
    find_common = squeeze(all(any(bsxfun(@eq, all_labels, permute(unique_val, [2 3 1])),2),1));
    common_val = unique_val(find_common);

    %set any cells that are not present in all frames to 0 (background)
    for idx = 1:180
        cell_matrix = importdata([matrix_storage_folder,filesep,'filtered_cell_matrix_',shot_cat{index},'_',num2str(idx),'.mat']);
        pixelsList = regionprops(cell_matrix, 'PixelIdxList');
        new_cell_mat = cell_matrix;
        cell_labels = unique(cell_matrix);
        cell_labels = nonzeros(cell_labels);
        containment = ismember(cell_labels,common_val);
        for j = 1: length(containment)
            disp(j)
            if ~containment(j)
                indices = pixelsList(cell_labels(j)).PixelIdxList;
                [y,x] = ind2sub(size(cell_matrix), indices);
                correct_idx = sub2ind(size(cell_matrix),y,x);
                new_cell_mat(correct_idx)=0;
            end
        end
        save([matrix_storage_folder,filesep,'consistent_cell_matrix_',shot_cat{index},'_',num2str(idx),'.mat'],'new_cell_mat');
        fname = ['composite_aligned_' num2str(idx, '%.1d') '.tif'];   
        disp(fname)
        cmap_mat = colormap(jet(max(unique(new_cell_mat))));
        img = imread(fullfile(base_dir,shot_cat{index},image_subdir,fname));
        label_img = labeloverlay(img,new_cell_mat,'Colormap',cmap_mat); %overlay filtered cell label matrix w image
        filtered_label_bounds = boundarymask(new_cell_mat);
        result_img = labeloverlay(label_img,filtered_label_bounds,'Transparency',0,'Colormap',sky);
        imshow(result_img)
        imwrite(result_img,[label_im_dir, filesep,shot_cat{index},'_consistent_w_boundaries_',num2str(idx),'.tif'],"tif");

    end

    toc
    %store cell numbers that remain in mat for each frame
    filtered_label_dir = [label_dir,filesep,'data',filesep,shot_cat{index},filesep,'labels'];
    save([filtered_label_dir,filesep,'filtered_labels_consistent_',shot_cat{index},'.mat'],'common_val');

end




     