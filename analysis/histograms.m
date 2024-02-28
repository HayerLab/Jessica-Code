clear;clc;close all; 
base_folder = 'Z:\Jessica\tracking_code\label_matrix\single_nuclei\data';
matrix_folder = 'matrix';
add_path;

%the shots to be tracked, change this if diff shots are used
shots_row = {'3'}; shots_col = {'8'}; shots_sites = {'1','2'};
delimiter = {'_'};
%format shots number
shot_cat = strcat(shots_row,delimiter,shots_col,delimiter,shots_sites);

for index = 1:numel(shot_cat)
    labels = zeros([180,10000]);
    for i = 1:180
        fname = ['consistent_cell_matrix_', shot_cat{index},'_', num2str(i), '.mat'];
        matrix = importdata([base_folder, filesep, shot_cat{index}, filesep, matrix_folder, filesep,fname]);
        values = nonzeros(unique(matrix))';
        labels(i,1:length(values)) = values;
        
    end
    labels(:,all(labels == 0))= [];
    unique_labels = nonzeros(unique(labels));
    occurences = arrayfun(@(x) sum(labels(:) == x), unique_labels);
    figure; h =histogram(occurences,36);
    xlabel({'Track length', '(frames)'});
    ylabel({'Frequency'});
    %h.BinLimits=[1 180];
end