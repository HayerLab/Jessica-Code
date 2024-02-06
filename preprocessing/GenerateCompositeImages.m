%%This Matlab scripts combine the DAPI and YFP microscopy images to create
%%an 8 bit RGB image with YFP as green channel and DAPI as red
%%channel.
%%Output: composite/merged channel image in .tif format for both the raw
%%images and the aligned raw images. 
%Matlab File Dependencies: Alignmentparameter.m;
%%%%%%%%%GenerateReferenceBackgroundImages.m;export_fig.m;
%%%%%%%%%calculate_bg_img_rm_blobs.m; showImagesMergeChannels.m;
%%%%%%%%%%%%%%%dualviewAlignFromFittedSurface.m

clear;clc;close all; set(0,'DefaultFigureVisible','off');

%set directories for root, raw images, backgriund, merged images
root = 'Z:\Jessica';
rawdir=[root,filesep,'raw',filesep,'160519_CDH5_mCitrine_Hoechst_1min_interval_20x',filesep, '3_7_1'];
bgdir=[root,filesep,'background'];
merged_channel_dir = [root, filesep,'raw', filesep, 'composites', filesep, '3_7_1'];
if ~exist(merged_channel_dir, 'dir')
    mkdir(merged_channel_dir)
end
load([bgdir,filesep,'3_7_1', filesep, 'alignment parameters pX pY.mat']);
position='3_7_1';
warning('off')

%% Load and process background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings

%background images for raw data
DAPIbg_raw=double(imread([bgdir,filesep,'AVG_rawdata_DAPI.tif']));
YFPbg_raw=double(imread([bgdir,filesep,'AVG_rawdata_YFP.tif']));
bg1(:,:,1)=DAPIbg_raw; 
bg1(:,:,2)=YFPbg_raw;

%%aligning bg images of the two channels
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
DAPIbg=bg2(:,:,1);
YFPbg=bg2(:,:,2);

for i = 1:60
    imRatio_raw={};maskFinal={};cellCoors={};
    YFPmaskFinal = {}; YFPcellCoors={};
    frameNum=i;
    disp([position,'__',num2str(frameNum)]);
    imDAPI_raw=double(imread([rawdir,filesep,position,'_DAPI_',num2str(frameNum),'.tif']));
    imYFP_raw=double(imread([rawdir,filesep,position,'_YFP_',num2str(frameNum),'.tif']));

    %%%%%% Align DAPI/YFP images
    imstack(:,:,1)=imDAPI_raw; imstack(:,:,2)=imYFP_raw;
    imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
    imDAPI_al=imaligned(:,:,1);
    imYFP_al=imaligned(:,:,2);
    i = i+1;

%% Display alignment
alignmentFigure = figure;
merged = subplot(1,2,1);
composite_file = [merged_channel_dir,filesep,'composite_',num2str(frameNum),'.tif'];
composite_file_aligned = [merged_channel_dir,filesep,'composite_aligned_',num2str(frameNum),'.tif'];
composite_dir = [merged_channel_dir,filesep,'raw_composites'];
composite_dir_aligned = [merged_channel_dir,filesep,'aligned'];

%create directories to store RGB images for raw images and aligned images
%if not already exist
mkdir ..\segmentation\3_7_1\Composites
mkdir ..\segmentation\3_7_1\Composites\raw_composites
mkdir ..\segmentation\3_7_1\Composites\aligned

%generates composite image with DAPI(nucleus) image as red channel and YFP
%(CDH5-mCitrine) as green channel
showImagesMergeChannels(composite_file,imstack(:,:,1),imstack(:,:,2)); %convert to 8-bit RGB images
newFig = figure;
tif_merged = copyobj(merged, newFig);
set(tif_merged, 'Position', get(0, 'DefaultAxesPosition'));
%print(gcf, '-dtiff',composite_file);
export_fig composite_file;

merged_al = subplot(1,2,2);
showImagesMergeChannels(composite_file_aligned,imDAPI_al,imYFP_al);
newFig_al = figure;
tif_merged_al = copyobj(merged_al, newFig_al);
set(tif_merged_al, 'Position', get(0, 'DefaultAxesPosition'));
%print(gcf, '-dtiff',composite_file_aligned);
export_fig composite_file_aligned;
end
