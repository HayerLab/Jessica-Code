function cell_num_track = cell_number(label_mat, frame_num, cell_label, track_mat, save)
% This function relabels the orignal cellpose label matrix with the
% corresponding cell number for each label based on tracking results.
% cell label matrix is cell matrix of size label_mat(1) x label_mat(2) if
%   saving each relabeled_matrix (for downstream filtering). 
% cell label matrix is matrix of size label_mat(1) x label_mat(2) x 180
%   if save = 0 and function called for use in approx. tracking accuracy

% %convert to cell matrix 
% cell_label_mat = mat2cell(resized_label_mat,size(label_mat,1), size(label_mat,2), 2);
cell_num = track_mat(:,frame_num,13); %cell number from mask tracking
cellID = 1:size(cell_num,1);
cellID =find(~isnan(cell_num));
mask_cell_match = cell_num(~isnan(cell_num));
tic
if save == 1 %only for saving each cell number label matrix
    for i = 1:length(mask_cell_match)
            [row,col] = find(label_mat==mask_cell_match(i)); %get indices of masks corresponding to cell
            label_idx = sub2ind([size(cell_label,1) size(cell_label,2)],row,col);
            % label cell matrix with cell number
            if isnan(cell_label(label_idx)) & cell_label(label_idx)~= 0
               ID = cellID(i);
               cell_label(label_idx) = ID;
            elseif ~isnan(cell_label(label_idx))
                cell_label(label_idx) = 0;
            elseif cell_label(label_idx) == 0
                continue;
            end
     end

else
    for i = 1:length(mask_cell_match)
        [row,col] = find(label_mat==mask_cell_match(i));
        label_idx = sub2ind([size(cell_label,1) size(cell_label,2) size(cell_label,3)], ...
            row,col,zeros(size(row,1),size(row,2))+frame_num);
        if isnan(cell_label(label_idx))
           ID = cellID(i);
           disp(ID);
           cell_label(label_idx) = ID;
           
        end
    end
end
toc
clear mask_cell_match
clear label_idx
cell_num_track = cell_label;

