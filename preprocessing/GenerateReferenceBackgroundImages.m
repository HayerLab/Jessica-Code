clear; clc;

% Specify folder where raw background images are located:
rawpath='C:\Users\Jessica Zhu\Downloads\raw\160519_CDH5_mCitrine_Hoechst_1min_interval_20x\3_8_2'; 
background='C:\Users\Jessica Zhu\Downloads\background\3_8_2';
% Specify the channels that were used, and which are part of the filenames
channels={'DAPI' 'YFP'};

for chan=1:numel(channels);
    calculate_bg_img_rm_blobs(rawpath,background,channels{chan});
    disp(channels(chan));
    disp(num2str(chan));
end

disp('done!');