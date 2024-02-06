function productcheck_getnucmask()
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_20x_Steve\';  %folder containing the movies folders [CHANGE]
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
nucr=8; %MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
%%%% sample specification %%%%%%%%%%%
row=1;
col=11;
site=1;
track=355;
%%% setup tempwell %%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
load([datadir,shot,'_alldata_FCT'],'bestsp','best_rc','corruptlist','leaveoutlist');
%%%%%%%%%%%%%%%%%%%%%%
crop = 1;   %0: no cropping, 1: crop
windowsize = 100;
%windowsize = 30;
%windowsize = 10;
windowhalf=windowsize/2;

%%% set up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempread=imread([cpdir,shot,'_nucedge&ring_',num2str(1),'.tif']);
[height,width]=size(tempread);
blocknum=1;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for f=1
    time1=tic;
    %%% read in pre-cropped images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_or=single(imread([cpdir,shot,'_nuc_',num2str(f),'.tif']));
    REs_or=single(imread([cpdir,shot,'_hDHB_',num2str(f),'.tif']));
    %%%%%  focus on track   %%%%%%%%%%%%%%
    if crop ==1
        cx=int16(bestsp{f}(track,2)); cx=cx-1*nucr;
        cy=int16(bestsp{f}(track,1));
        miny=cy-(windowhalf-1)*nucr; maxy=cy+windowhalf*nucr;
        minx=cx-(windowhalf-1)*nucr; maxx=cx+windowhalf*nucr;
        if minx<1
            minx=1; maxx=1+windowsize*nucr;
        end
        if miny<1
            miny=1; maxy=1+windowsize*nucr;
        end
        if maxx>height
            maxx=height; minx=height-windowsize*nucr;
        end
        if maxy>width
            maxy=width; miny=width-windowsize*nucr;
        end
        DAs_or=DAs_or(minx:maxx,miny:maxy);
        REs_or=REs_or(minx:maxx,miny:maxy);
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_bs=bgsub_MC(DAs_or,blockheight,blockwidth);
    
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_pad=product_getnucmask_rev01(DAs_bs,nucr);
    
    %%% get nuclear edges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nucmoviering=bwmorph(DAs_pad,'remove');

    %%% build image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tempframe=imadjust(mat2gray(DAs_bs));
    tempframe(:,:,2)=imadjust(mat2gray(nucmoviering));
    %tempframe(:,:,3)=imadjust(mat2gray(nucmoviering));
    tempframe(:,:,3)=0;

    imshow(imresize(tempframe,2));

    toc(time1)
end

cd ..\Analysis; %return to this directory

end