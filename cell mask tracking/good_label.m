function good_labels = good_label(unique_labels,cell_mat,size_tuple)
%function to determine which labels are at the edge of the image matrix
%good_labels is a struct array with the x, y coordinates of labels not at
%edge

pixels = regionprops(cell_mat, 'PixelIdxList');
structure_size = length(pixels);
good_labels(structure_size) = struct();

for i = 1:length(pixels)
    [y,x] = ind2sub(size(cell_mat), pixels(i).PixelIdxList);
    checkEdgeY = any(y(:) <= 5 |y(:) ==  size(cell_mat,1));
    checkEdgeX = any(x(:) <= 5 |x(:) ==  size(cell_mat,2));
    if checkEdgeY ~= true && checkEdgeX ~= true
        good_labels(i).x = x;
        good_labels(i).y = y;
    end
end
