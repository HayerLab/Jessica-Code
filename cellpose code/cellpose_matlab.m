%%%%
% This is the MATLAB script for batch segmented microscopy images
% it is the equivalent of the python script

%%%%%%%% Set up folders
drive = "Z:\";
user = "Jessica\";
base_folder = "segment\";
data_folder = "matrix\"; %storage location of cellpose label matrix
model_filepath = "model with varied sets\training_set\models\";
image_folder = "Composites\aligned\tif\"; %folder containing RGB image ot be segmented
img_folder = "images"; %for storing overlaid images
shots_row = {'3'}; shots_col = {'8'}; shots_sites = {'1','2'};
delimiter = {'_'};
%format shots number
shot_cat = strcat(shots_row,delimiter,shots_col,delimiter,shots_sites);

%Specify model path if using custom trained model
model_path = strcat(drive,user,model_filepath);
diameter = 104;
cp = cellpose(Model="HUVEC_monolayer_varied_set",ModelFolder=model_path);

for index = 1:numel(shot_cat)
    file_path = strcat(drive,user,base_folder,shot_cat{index},filesep,image_folder);
    files = dir(fullfile(file_path, '*.tif'));
    files = natsortfiles(files); %sort files in natural order
    files = struct2cell(files);
    files = files(1,:);
    for i = 1:length(files)
        img = imread(strcat(drive,user,base_folder,shot_cat{index},filesep,image_folder,files(i)));
        % G channel use as main channel, R channal as auxiliary (nuclear)
        masks = segmentCells2D(cp,img(:,:,2),ImageCellDiameter=diameter,AuxiliaryChannelImage=img(:,:,1));
        save(strcat(drive,user,base_folder,shot_cat{index},filesep,data_folder,filesep,num2str(i),".mat"),'masks');
        labels = labeloverlay(img,masks);
        imshow(labels)
        imwrite(img,strcat(drive,user,base_folder,shot_cat{index},filesep,img_folder,filesep,num2str(i),'.tif'),"tif");

    end
    
end
