clc
clear
close all

label_dir = 'Z:\Jessica\tracking_code\label_matrix';
label_im_dir = 'Z:\Jessica\tracking_code\label_matrix\images';
base_dir = 'Z:\Jessica\segment'; %base directory for segmented masks
image_subdir = '\Composites\aligned\tif\';
shots_row = {'3'}; shots_col = {'8'}; shots_sites = {'1','2'};
delimiter = {'_'};
%format shots number
shot_cat = strcat(shots_row,delimiter,shots_col,delimiter,shots_sites);


tic

for index = 1:numel(shot_cat)
    shot_num = shot_cat{index};
    disp(shot_num)

    if index > 1
        clear filtered_labels_list
    end 

    bound_im_dir = [label_im_dir,filesep,'boundary_outline',filesep,shot_cat{index}];
    if ~exist(bound_im_dir,'dir')
        mkdir(bound_im_dir);
    end

    cell_im_dir = [label_im_dir,filesep,'cell_labels',filesep,shot_cat{index}];
    if ~exist(cell_im_dir,'dir')
        mkdir(cell_im_dir);
    end

    matrix_storage_folder = [label_dir,filesep,'data',filesep,shot_cat{index},filesep,'matrix'];
    if ~exist(matrix_storage_folder,'dir')
        mkdir(matrix_storage_folder);
    end

    filtered_label_dir = [label_dir,filesep,'data',filesep,shot_cat{index},filesep,'labels'];
    if ~exist(filtered_label_dir,'dir')
        mkdir(filtered_label_dir);
    end

    label_image_dir = [label_im_dir,filesep,'filtered',filesep,shot_cat{index}];
    if ~exist(label_image_dir,'dir')
        mkdir(label_image_dir);
    end

    filtered_labels_list(180) = struct();
    for k = 1:180
        cell_mat = importdata([label_dir,filesep,shot_cat{index},filesep,strcat(shot_cat{index},'_',num2str(k),'_cell_num_matrix.mat')]);
        labels = unique(cell_mat);
        label_bounds = boundarymask(cell_mat);
        filtered_labels = good_label(labels,cell_mat,size(cell_mat));
        filtered_cell_num = find(~cellfun(@isempty,{filtered_labels.x}));
        filtered_labels_list(k).cell_label = filtered_cell_num;
        fname = ['composite_aligned_' num2str(k, '%.1d') '.tif'];   
        disp(fname)
        img = imread(fullfile(base_dir,shot_cat{index},image_subdir,fname));
        cmap = colormap(jet(size(filtered_cell_num,2)));
        boundary_mask = labeloverlay(img,label_bounds,'Transparency',0);
       
        imwrite(boundary_mask,[bound_im_dir, filesep,shot_cat{index}, '_boundary_',num2str(k),'.tif'],"tif");
        %imshow(boundary_mask)

        cmap_mat = colormap(jet(size(filtered_labels,2)));
        cell_overlay = labeloverlay(img,cell_mat,'Colormap',cmap_mat);
        %imshow(cell_overlay);
        
        imwrite(cell_overlay,[cell_im_dir, filesep,shot_cat{index}, '_cell_labels_',num2str(k),'.tif'],"tif");
        filtered_label_im = labeloverlay(img,cell_mat,'Colormap',cmap,'IncludedLabels',filtered_cell_num);
        %figure; imshow(filtered_label_im)
        pos=1:length(filtered_labels);
        pos(filtered_cell_num)=[];
        pixelsList = regionprops(cell_mat, 'PixelIdxList');
        new_cell_mat = cell_mat;
        for j = 1: length(pos)
            disp(j)
            indices = pixelsList(pos(j)).PixelIdxList;
            [y,x] = ind2sub(size(cell_mat), indices);
            correct_idx = sub2ind(size(cell_mat),y,x);
            new_cell_mat(correct_idx)=0;
        end
        
        save([matrix_storage_folder,filesep,'filtered_cell_matrix_',shot_cat{index},'_',num2str(k),'.mat'],'new_cell_mat');
        %imshow(filtered_label_im)
        %figure;
        %imshow(labeloverlay(img,new_cell_mat))
        
        imwrite(filtered_label_im,[label_im_dir, filesep,shot_cat{index},'_filtered_',num2str(k),'.tif'],"tif");
    end 
    
    
    save([filtered_label_dir,filesep,'filtered_labels_',shot_cat{index},'.mat'],'filtered_labels_list');
end
toc


