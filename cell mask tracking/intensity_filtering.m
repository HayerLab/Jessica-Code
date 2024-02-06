clc; clear; close all;
%%%% This script is the second step for cell mask track quality control.
%%%% Run this script after running remove_border.m
%%%% This script filters out any cells whose mean intensity exceeds the
%%%%    average mean intensity (avg across all frames) beyond + one std

%%%%%% set up paths %%%%%%%
label_dir = 'Z:\Jessica\tracking_code\label_matrix';
label_im_dir = 'Z:\Jessica\tracking_code\label_matrix\images';
base_dir = 'Z:\Jessica\segment'; %base directory for segmented masks
image_subdir = '\Composites\aligned\tif\';
rawdir = 'Z:\Jessica\raw\160519_CDH5_mCitrine_Hoechst_1min_interval_20x';
shots_row = {'3'}; shots_col = {'8'}; shots_sites = {'1','2'};
delimiter = {'_'};
%format shots number
shot_cat = strcat(shots_row,delimiter,shots_col,delimiter,shots_sites);

for index = 1%2:numel(shot_cat)
    %path where cell matrix is stored
    matrix_storage_folder = [label_dir,filesep,'data',filesep,shot_cat{index},filesep,'matrix'];
    if index > 1
        clear all_means
    end 
    all_means(180) = struct();

    %loop through all frames and store intensity properties of labels in
    %each frame in a struct array
    for i = 1:180
        disp(i)
        %import the cell matrix output with cells touching image border removed
        cell_matrix = importdata([matrix_storage_folder,filesep,'filtered_cell_matrix_',shot_cat{index},'_',num2str(i),'.mat']);
        %raw_YFP = imread([rawdir,filesep,shot_cat{index},filesep,shot_cat{index},'_YFP_',num2str(i),'.tif']);
        
        %read composite image and convert to grayscale image
        composite = imread([base_dir,filesep,shot_cat{index},filesep,image_subdir,filesep,'composite_aligned_',num2str(i),'.tif']);
        grayscale = rgb2gray(composite);
        
        %get intensity properties using reigonprops
        labels_intensity = regionprops('table', cell_matrix,grayscale,'MeanIntensity','MaxIntensity','MinIntensity');
        meanI = labels_intensity.MeanIntensity;
        padded = NaN([1000 1]); %NaN array so array of intensities same length across all frames
        padded(1:length(meanI)) = meanI;
        maxI = labels_intensity.MaxIntensity;
        minI = labels_intensity.MinIntensity;
        all_means(i).mean = padded;
        % all_means(i).maxIntensity = maxI;
        % all_means(i).minIntensity = minI;

    end
    
    all_params = cell2mat(struct2cell(all_means));
    avg_mean = mean(all_params,3,'omitmissing'); %avg of mean intensity for each label across all frames
    threshold = mean(all_params,'all','omitmissing')+std(all_params,[],'all','omitmissing');
    thresholded_cell = find(avg_mean> threshold); %cells exceeding threshold to be removed
    filtered_labels_list(180) = struct();
    for idx = 1:180
        %import cell matrix data
        cell_matrix = importdata([matrix_storage_folder,filesep,'filtered_cell_matrix_',shot_cat{index},'_',num2str(idx),'.mat']);
        composite = imread([base_dir,filesep,shot_cat{index},filesep,image_subdir,filesep,'composite_aligned_',num2str(idx),'.tif']);
        new_cell_mat = cell_matrix;

        %get idx of all pixels in each cell label in matrix
        pixelsList = regionprops(cell_matrix,'PixelIdxList'); 
        cell_num = unique(cell_matrix);
        %for loop for removing cells with high intensity by setting to 0
        %(background)
        for j = 1: length(thresholded_cell)
            %disp(j)
            if any(cell_num==thresholded_cell(j))
                indices = pixelsList(thresholded_cell(j)).PixelIdxList;
                [y,x] = ind2sub(size(cell_matrix), indices);
                correct_idx = sub2ind(size(cell_matrix),y,x);
                new_cell_mat(correct_idx)=0;
             end
           
        end
        filtered_cell_num = unique(new_cell_mat(new_cell_mat~=0));
        filtered_labels_list(idx).cell_label = filtered_cell_num;
        img = labeloverlay(composite,new_cell_mat); %overlay filtered cell label matrix w image

        %store overlaid image and label matrix data
        imwrite(img,[label_dir,filesep,'data',filesep,shot_cat{index},filesep,'no_clusters',filesep,'images',filesep,shot_cat{index},'_final_filtered_',num2str(idx),'.tif'],"tif");
        %figure;imshow(img)
        save([label_dir,filesep,'data',filesep,shot_cat{index},filesep,'no_clusters',filesep,shot_cat{index},'_final_filtered_',num2str(idx),'.mat'],'new_cell_mat');
    end
    
    %store cell numbers that remain in mat for each frame
    filtered_label_dir = [label_dir,filesep,'data',filesep,shot_cat{index},filesep,'labels'];
    save([filtered_label_dir,filesep,'filtered_labels_no_clusters_',shot_cat{index},'.mat'],'filtered_labels_list');
end


