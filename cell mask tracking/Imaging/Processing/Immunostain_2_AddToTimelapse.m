function Immunostain_2_AddToTimelapse(row,col,site)
%row=2;col=7;site=1;
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='H:\Documents\Projects\';
imagepath='E:\';
experimentpath='20130607 p21dCy2\20130719 EdU CycD1 Drug&siRNA panel\';
shadingpath='H:\Images\ShadingImages\20140424 DCYTC 10x\';

shot=wellnum2str(row,col,site);
%shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%datadir=([projectpath,experimentpath,'Data\']);
datadir=([projectpath,experimentpath,'Data_20140606\']);
separatedirectories=1;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask',shot,'_'];
end
maskwrite=0;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name1='CFP_';
nucedgename='nucedge_';
name2='YFP_';
name3='TexasRed_';
IFname1='DAPI_';
IFname2='Cy5_';
moviebin=1;
if moviebin==1
    nucr=12;
    debrisarea=100; %MCF-10A: 100
    %nucr=20;
    %debrisarea=600; %MCF-10A 20xBin1
elseif moviebin==2
    nucr=6;
    debrisarea=50; %MCF-10A:100, BJ5:50
end
boulderarea=20*debrisarea;
%blobthreshold=-0.02; %MCF10A 10xBin1
blobthreshold=-0.03; %MCF10A 20xBin1
timetotal=tic;
%%% load timelapse data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
[totalcells,totalframes,totalsignals]=size(tracedata);
%%% get previous mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawprev=single(imread([rawdir,name1,num2str(totalframes),'.tif']));
[nuc_mask_prev,~]=blobdetector_foreground(log(rawprev),nucr,blobthreshold,debrisarea);
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([shadingpath,'BG','.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
load([imagepath,experimentpath,'Raw\IBG\CFP_stain_',num2str(site),'.mat']); sc1=inferredbg;
load([imagepath,experimentpath,'Raw\IBG\YFP_stain_',num2str(site),'.mat']); sc2=inferredbg;
load([imagepath,experimentpath,'Raw\IBG\TexasRed_stain_',num2str(site),'.mat']); sc3=inferredbg;
load([imagepath,experimentpath,'Raw\IBG\DAPI_stain_',num2str(site),'.mat']); sc4=inferredbg;
load([imagepath,experimentpath,'Raw\IBG\Cy5_stain_',num2str(site),'.mat']); sc5=inferredbg;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw1=single(imread([rawdir,name1,'stain.tif'])); raw1=(raw1-bgcmos)./sc1;
raw2=single(imread([rawdir,name2,'stain.tif'])); raw2=(raw2-bgcmos)./sc2;
raw3=single(imread([rawdir,name3,'stain.tif'])); raw3=(raw3-bgcmos)./sc3;
IFraw1=single(imread([rawdir,IFname1,'stain.tif'])); IFraw1=(IFraw1-bgcmos)./sc4;
IFraw2=single(imread([rawdir,IFname2,'stain.tif'])); IFraw2=(IFraw2-bgcmos)./sc5;
NEfile=[maskdir,nucedgename,'stain.tif'];
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nuc_mask=blobdetector(log(raw1),nucr,blobthreshold,debrisarea);
%[nuc_mask,foreground]=blobdetector_foreground(log(raw1),nucr,blobthreshold,debrisarea);
nuc_mask=blobdetector_3(log(raw1),nucr,blobthreshold,debrisarea);
foreground=nuc_mask;
nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea);
%nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
%%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[height,width]=size(raw1);
nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
border=bwareaopen(nuc_mask,height*2+width*2-4);
nuc_mask=logical(nuc_mask-border);
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanmask=imdilate(foreground,strel('disk',nucr/2));
nanmaskcyto=imdilate(foreground,strel('disk',nucr));
blur1=imfilter(raw1,fspecial('disk',3),'symmetric');
real1=bgsubmasked_3(blur1,nanmask,1);
blur2=imfilter(raw2,fspecial('disk',3),'symmetric');
real2=bgsubmasked_3(blur2,nanmaskcyto,1);
blur3=imfilter(raw3,fspecial('disk',3),'symmetric');
real3=bgsubmasked_3(blur3,nanmask,1);
IFblur1=imfilter(IFraw1,fspecial('disk',3),'symmetric');
IFreal1=bgsubmasked_3(IFblur1,nanmask,1);
IFblur2=imfilter(IFraw2,fspecial('disk',3),'symmetric');
IFreal2=bgsubmasked_3(IFblur2,nanmask,1);
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nuc_label,numcells]=bwlabel(nuc_mask);
nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
%%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_mass=nuc_density.*nuc_area;
%%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reljitx,reljity]=detectGlobalShift(nuc_mask_prev,nuc_mask);
reljitx=-reljitx;
reljity=-reljity;
reljitter=[reljitx,reljity];
prevjitter=jitters(totalframes,:);
IFjitter=prevjitter+reljitter;
nuc_center(:,1)=nuc_center(:,1)+IFjitter(1);
nuc_center(:,2)=nuc_center(:,2)+IFjitter(2);
%%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
%%% track & correct merges (update centers, masses and labels) %%%%%%%%%%%%
debugpackage={rawprev,nuc_mask_prev,nuc_mask,prevjitter,reljitter};
[tracked,nuc_label]=adaptivetrack_IF(tracedata(:,totalframes,1:4),initdata,nuc_label,nucr,debugpackage);
numcells=size(tracked,1);
nuc_center=tracked(:,[1 2]);
nuc_area=tracked(:,3);
nuc_mass=tracked(:,4);
%%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if maskwrite
    extractmask=bwmorph(nuc_label,'remove');
    imwrite(uint16(extractmask),NEfile);
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellid=find(~isnan(nuc_mass));
nuc_info=regionprops(nuc_label,'PixelIdxList');
ring_label=getcytoring_3(nuc_label,4,real2); %10xB1:4 20xB1:4
ring_info=regionprops(ring_label,'PixelIdxList');
nanvec=ones(numcells,1)*NaN; sig1=nanvec; sig2=nanvec; sig3=nanvec;
sig2ring_75th=nanvec; sig2ring_mode=nanvec;
IFsig1=nanvec; IFsig2=nanvec;
for cc=cellid'
    sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
    sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
    sig3(cc)=median(real3(nuc_info(cc).PixelIdxList));
    IFsig1(cc)=mean(IFreal1(nuc_info(cc).PixelIdxList));
    IFsig2(cc)=mean(IFreal2(nuc_info(cc).PixelIdxList));
    
    ringall=real2(ring_info(cc).PixelIdxList); 
    ringall(ringall>prctile(ringall,95))=[]; %hist(ringall,25);
    sig2ring_75th(cc)=prctile(ringall,75); %previously mean

    truering=ringall(ringall>25);

    if numel(truering)<50
         truering=ringall;
    end
    if numel(ringall)<2
         sig2ring_mode(cc)=NaN;
         continue;
    end             

    bmin=min(truering); bmax=max(truering); bstep=(bmax-bmin)/25; bins=bmin:bstep:bmax;
    [kval,xval]=ksdensity(truering,bins);
    maxidx=find(kval==max(kval),1); %first mode
    sig2ring_mode(cc)=xval(maxidx);
end
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([imagepath,experimentpath,'Raw\','bleedthroughrate_RFPtoCy5.mat'],'bleedthroughrate');
IF_bleedthrough=bleedthroughrate(1)+sig3*bleedthroughrate(2);
IFsig2=IFsig2-IF_bleedthrough;
%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig2,sig3,sig2ring,IFsig1,IFsig1ring,IFsig2];
%IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,sig2ring,IFsig1];
IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,IFsig1,IFsig2,sig2ring_75th,sig2ring_mode];
save([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
toc(timetotal);
%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(IFraw1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
imshow(tempframe);
%}