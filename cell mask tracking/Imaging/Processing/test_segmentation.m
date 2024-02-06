%% test segmentation
clear;clc;
rawdir='K:\Rachel\161103_IXM\raw\';

nucr=15;
blobthreshold=-0.02;
debrisarea=350;

raw1=single(imread([rawdir,'4_6_2_DAPI_50.tif']));

nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);

% merge=zeros(size(raw1,1),size(raw1,2),3);
% merge(:,:,1)=nuc_mask;
% merge(:,:,2)=mat2gray(raw1,[600 8000]);
% 
% imshow(merge);


extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(raw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure(nucr),imshow(tempframe);


%% 
% clear;
% nucr=3;
% rawdir='D:\data\140708_Axon\Mitofission\';
% raw1=single(imread([rawdir,'2_1_1_DAPI_1.tif']));
% path(path,'D:\Stanford\Matlab\Tracking_Min_140611\Imaging\Functions\General');
% DAs_bl=imfilter(log(raw1),fspecial('disk',round(nucr/1)),'symmetric');
% DAs_bs=bgsub(DAs_bl,3*nucr,0);
% nuc_mask=getdapimask(DAs_bs,nucr);
% extractmask=bwmorph(nuc_mask,'remove');
% tempframe=imadjust(mat2gray(raw1));
% tempframe(:,:,2)=extractmask;
% tempframe(:,:,3)=0;
% figure,imshow(tempframe);