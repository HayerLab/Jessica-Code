row=3;col=4;site=1;track=63;
%shot=wellnum2str(row,col,site);
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
%%% set filepaths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imagepath='E:\';
imagepath='H:\Images\';
projectpath='H:\Documents\Projects\';
%experimentpath='20130607 p21dCy2\20130719 EdU CycD1 Drug&siRNA panel\';
experimentpath='Michael\';

rawdir=[imagepath,experimentpath,'BGsubbed\',shot,'\'];
maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
datadir=([projectpath,experimentpath,'Data\']);
namemask='nucedge_';
name1='CFP_real_';
name2='YFP_real_';
%name3='TexasRed_real_';
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy);
%%% set viewing parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
winrad=50;
imagescale=5;

shottrack=[shot,'_track',num2str(track)];
M=VideoWriter([datadir,shottrack],'Uncompressed AVI');
M.FrameRate = 4;
open(M);
tempread=imread([maskdir,namemask,num2str(1),'.tif']);
[height,width]=size(tempread);

for f=tracestats(track,1):tracestats(track,2)
    %%% read in pre-cropped images %%%%%%%
    nucmask=single(imread([maskdir,namemask,num2str(f),'.tif']));
    raw1=single(imread([rawdir,name1,num2str(f),'.tif']));
    raw2=single(imread([rawdir,name2,num2str(f),'.tif']));
    %raw3=single(imread([rawdir,name3,num2str(f),'.tif']));
    %%%%%  focus on track   %%%%%%%%%%%%%%
    cx=int16(tracedata(track,f,1)-jitters(f,1));
    cy=int16(tracedata(track,f,2)-jitters(f,2));
    miny=cy-winrad; maxy=cy+winrad;
    minx=cx-winrad; maxx=cx+winrad;
    if minx<1
        minx=1; maxx=1+2*winrad;
    end
    if miny<1
        miny=1; maxy=1+2*winrad;
    end
    if maxx>width
        maxx=width; minx=width-2*winrad;
    end
    if maxy>height
        maxy=height; miny=height-2*winrad;
    end
    nucmask=nucmask(miny:maxy,minx:maxx);
    raw1=raw1(miny:maxy,minx:maxx);
    raw2=raw2(miny:maxy,minx:maxx);
    %raw3=raw3(miny:maxy,minx:maxx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Add frame to movie
    tempframe=mat2gray(imresize(raw1,imagescale));
    tempframe(:,:,2)=mat2gray(imresize(raw2,imagescale));
    tempframe(:,:,3)=0;
    tempframe(:,:,1)=mat2gray(imresize(nucmask,imagescale));
    
    writeVideo(M,im2frame(tempframe));
end
close(M);
