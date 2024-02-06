function mask_tracking = mask_matching(file_path, tracedata,shot)
files = dir(fullfile(file_path, '*.csv'));
files = natsortfiles(files); %sort files in natural order (alphanumeric)
files = struct2cell(files);
files = files(1,:);

%resize tracedata matrix
new_size = [size(tracedata,1), size(tracedata,2), size(tracedata,3)+7];
resized_trace = NaN(new_size);
resized_trace(1:size(tracedata,1), 1:size(tracedata,2), 1:size(tracedata,3)) = tracedata;
x_coord = round(tracedata(:,:,1),0); %x and y coordinate of nucleus centroid from nuclear tracking
y_coord = round(tracedata(:,:,2),0);

for i = 1:length(files)
    disp(i)
    tic
    %read mask csv file to retrieve label matrix for cell masks
    mask = readmatrix(strcat(file_path, filesep, files{i}));

    %calculate cell shape params from masks
    mask_info = struct2cell(regionprops(mask,'Area','Centroid', 'Circularity', 'Eccentricity','Perimeter')');
    mask_area = squeeze(cell2mat(mask_info(1,1,:)));
    mask_centroid = squeeze(cell2mat(mask_info(2,1,:)'));
    mask_circ = squeeze(cell2mat(mask_info(3,1,:)));
    mask_ecc = squeeze(cell2mat(mask_info(4,1,:)));
    mask_perim = squeeze(cell2mat(mask_info(5,1,:)));

    %match coords of nuclei centre to mask
    indices = sub2ind([size(mask,1) size(mask,2)], y_coord(:,1),x_coord(:,1)); %convert to lin indexing
    nucleiID = 1:size(x_coord,1); %vector to keep track of nuclei IDs
    nucleiID = nucleiID(~isnan(indices)); %remove nuclei IDs with NaN value
    
    values = mask(indices(nucleiID)); %Mask value at each nucleus position
    nucleiID = nucleiID(values~=0);
    values = mask(indices(nucleiID));

    %fill added matrix layers with mask tracking info & params
    resized_trace(nucleiID, i, size(tracedata,3)+1) = values;
    resized_trace(nucleiID, i, size(tracedata,3)+2) = mask_area(values);
    resized_trace(nucleiID, i, size(tracedata,3)+3) = mask_centroid(values,1);
    resized_trace(nucleiID, i, size(tracedata,3)+4) = mask_centroid(values,2);
    resized_trace(nucleiID, i, size(tracedata,3)+5) = mask_circ(values);
    resized_trace(nucleiID, i, size(tracedata,3)+6) = mask_ecc(values);
    resized_trace(nucleiID, i, size(tracedata,3)+7) = mask_perim(values);
    toc
end
                            
resized_trace(resized_trace==0) = NaN;
mask_tracking = resized_trace; %output new tracking matrix
save([file_path,'mask_tracking_',shot,'.mat'],'resized_trace'); %save matrix
