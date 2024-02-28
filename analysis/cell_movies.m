nImage = 180;                                                        % L1
fps = 10.0;                                                           % L2
%image_folder = 'Z:\Jessica\tracking_code\label_matrix\single_nuclei\images\filtered'; %set file location for images
%image_folder = '/Volumes/hayerlab2/Jessica/segment/3_8_2/Composites/aligned/tif';
image_folder = 'Z:\Jessica\tracking_code\label_matrix\single_nuclei\images\consistent\3_8_2';
video_folder = 'Z:\Jessica\tracking_code\label_matrix\single_nuclei\movie'; 
if ~exist(video_folder,'dir')
   mkdir(video_folder);
end
%video_folder = '/Volumes/hayerlab2/Jessica/segment/3_8_2/Composites/aligned/movie';
oVideo = VideoWriter(fullfile(video_folder, '3_8_2_consistent.mp4'));       % L5
oVideo.FrameRate = fps;                                              % L6
open(oVideo)   

rng('default')
rng(2)
num = randperm(400);
cmap = colormap(hsv(400));  
randomized_cmap = cmap(num,:);
for i = 1:nImage                                                     % L8
    fname = ['3_8_2_consistent_w_boundaries_', num2str(i, '%.1d'), '.tif'];   
    disp(fname)
    img = imread(fullfile(image_folder,fname));
    imshow(img)
    curImage = img;
    writeVideo(oVideo, curImage);   
end                                                                  % L12
close(oVideo)   

