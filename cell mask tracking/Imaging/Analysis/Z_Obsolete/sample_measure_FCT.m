function sample_measure_org()
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_20x_Steve\';  %folder containing the movies folders [CHANGE]
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
nucr=8; %MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
movie=1;
singlecell=0;
%%%% sample specification %%%%%%%%%%%
row=1;      %Geminin_10x_Steve: FCT: 2_1_1_779
col=11;
site=1;
track=416;
%%% setup tempwell %%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
shottrack=[shot,'_track',num2str(track)];
load([datadir,shot,'_alldata_FCT'],'bestsp','best_rc','corruptlist','leaveoutlist');
%%%%%%%%%%%%%%%%%%%%%%
colorcode = jet(128);
colorcode1 = colorcode(:,1);
colorcode2 = colorcode(:,2);
colorcode3 = colorcode(:,3);
windowsize = 30;
%windowsize = 10;
windowhalf = windowsize/2;
cropsize = 8;
crophalf = cropsize/2;
difhalf = windowhalf-crophalf;
framezero = 0;   %Geminin_Steve_10x: 63

if movie
    M=VideoWriter([datadir,shottrack],'Uncompressed AVI');
    M.FrameRate = 4;
    open(M);
end
tempread=imread([cpdir,shot,'_nucedge&ring_',num2str(1),'.tif']);
[m,n]=size(tempread);

%for f=best_rc(track,1):best_rc(track,3)
for f=framezero+1
    time1=tic;
    %%% read in pre-cropped images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_or=single(imread([cpdir,shot,'_nuc_',num2str(f),'.tif']));
    REs_or=single(imread([cpdir,shot,'_hDHB_',num2str(f),'.tif']));
    %CEs_or=single(imread([tifdir,'/',movieName,'_geminin_',num2str(f),'.tif']));
    %%%%%  focus on track   %%%%%%%%%%%%%%
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
    if maxx>m
        maxx=m; minx=m-windowsize*nucr;
    end
    if maxy>n
        maxy=n; miny=n-windowsize*nucr;
    end
    DAs_or=DAs_or(minx:maxx,miny:maxy);
    REs_or=REs_or(minx:maxx,miny:maxy);
    %CEs_or=CEs_or(minx:maxx,miny:maxy);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% image processing
    DAs_bl=log(imfilter(DAs_or,fspecial('disk',floor(nucr/2)),'symmetric')); 
    DAs_bs=bgsub(DAs_bl,10*nucr,0.05);
    %DAs_bs_sharp=bgsub_test(DAs_or,10*nucr,0.05);
    %REs_bl=(imfilter(REs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    REs_bs=bgsub(REs_or,10*nucr,0.05);
    %CEs_bl=(imfilter(CEs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    %CEs_bs=bgsub(CEs_bl,10*nucr,0.05);
    
    %% get data
    DAs_ma=getdapimask(DAs_bs,nucr);  %make the mask. this segmentation function is the engine of Feng Chiao's script
    nuclearedges=bwmorph(DAs_ma,'remove');
    %DAs_la=bwlabel(DAs_ma);  %labels the objects with numbers
    dapimaskdilated=imdilate(DAs_ma, strel('disk', 1, 8));  %dilate dapimask by 2 pixels; 8 means octagon   
    ring=imdilate(dapimaskdilated, strel('disk', 3, 8)) - dapimaskdilated;  %dilate by 2 pixels, and subtract to make a ring              
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Add frame to movie
    if movie
        
        tempframe=imadjust(mat2gray(nuclearedges));
        tempframe(:,:,2)=imadjust(mat2gray(ring));
        tempframe(:,:,3)=imadjust(mat2gray(REs_bs));
        %tempframe(:,:,3)=0;
        %}
        
        %{
        %DAs_bs=imfilter(DAs_bs,fspecial('disk',1),'symmetric');
        DAs_bs = imadjust(mat2gray(DAs_bs));
        normimage = DAs_bs/(max(max(DAs_bs)));
        codeimage = round(normimage*128);
        codeimage(codeimage==0) = 1;
        tempframe = colorcode1(codeimage);
        tempframe(:,:,2) = colorcode2(codeimage);
        tempframe(:,:,3) = colorcode3(codeimage);
        %}
        
        %tempframe = tempframe(1+difhalf*nucr:end-difhalf*nucr,1+difhalf*nucr:end-difhalf*nucr,:);
        imshow(imresize(tempframe,2));
        %writeVideo(M,im2frame(tempframe));
    end
    toc(time1)
end

%save([datadir,shottrack,'attributes.mat'],'mon','msn','sbn','ssn','aor','anr','sor','snr','areabignuc','areasmallnuc');
if movie
    close(M);
    set(gcf,'PaperPosition',[0 0 3 3]);
    saveas(gcf,'h:\Downloads\Fig.jpg');
end
cd ..\Analysis; %return to this directory

end