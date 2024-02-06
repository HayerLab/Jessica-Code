%function B_extractfeatures_MC(row,col,site)
row=2;col=4;site=1;
timetotal=tic;
cd ..\Functions;
path = 'h:\Documents\Timelapse\Timescape\Debugging\';
datadir = path;
cpdir = path;
SF=17;EF=17; %Sabrina20x:208 Steve20x:218 Steve10x:110 Steve&Sabrina:240
initF=SF;   %the first intended frame (only matters if I have to restart
nucr=25;
DAname = '_nuc_';
%DAname = '_mRFP1_';
DAedgename = '_nucedge&ring_';
%REname = '_hDHB_';
%CEname = '_cdt1_';  %Pilot:'cdt1' Steve20x:'geminin' Steve10x:'geminin' Steve&Sabrina:'cdt1'
debrisrad = nucr/5;
minnucarea=pi*(debrisrad)^2;
continuation=0; restartframe=76;

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
if continuation==1
    load([datadir,'wellsss_',shot,'_restart'],'wellsss');
    wellsssrestart=cell(1,1,EF-SF+1);
    wellsssrestart(1:restartframe-1)=wellsss(1:restartframe-1);
    wellsss=wellsssrestart;
    SF=restartframe;
else
    wellsss=cell(1,1,EF-SF+1);  %row, column, frame (for 96-well plate)  
end
tempread=imread([cpdir,shot,DAname,num2str(17),'.tif']);
[height,width]=size(tempread);
blocknum=3;
blockheight=ceil(height/blocknum);
blockwidth=ceil(width/blocknum);

for f=SF:EF
    fprintf('frame %0.0f\n',f);
    timeframe=tic;

    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAfile=[shot,DAname,num2str(f),'.tif'];
    NEfile=[shot,DAedgename,num2str(f),'.tif'];
    %REfile=[shot,REname,num2str(f),'.tif'];
    %CEfile=[shot,CEname,num2str(f),'.tif'];  
    
    DAs_or=single(imread([cpdir,DAfile]));
    %REs_or=single(imread([cpdir,REfile]));
    %CEs_or=single(imread([cpdir,CEfile]));
    %%% image processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_bs=bgsub_MC(log(DAs_or),blockheight,blockwidth);
    %REs_bs=bgsub_MC(log(REs_or),blockheight,blockwidth);
    %CEs_bs=bgsub_MC(log(CEs_or),blockheight,blockwidth);
    %%% nuclear segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_pad=getnucmask_histsweep(DAs_bs,nucr);  %MC histogram sweep & concave detector
    %{
    midrad=2;
    ringwidth=midrad+1;
    [height,width]=size(DAs_pad);
    zeroborder=ones(height,width);
    zeroborder(1,:)=0;zeroborder(end,:)=0;zeroborder(:,1)=0;zeroborder(:,end)=0;
    ringzeroborder=ones(height,width);
    ringzeroborder(1:ringwidth+1,:)=0;ringzeroborder(end-ringwidth:end,:)=0;ringzeroborder(:,1:ringwidth+1)=0;ringzeroborder(:,end-ringwidth:end)=0;
    DAs_pad=DAs_pad & zeroborder;   %necessary to imerode from edges
    DAs_pad=imopen(DAs_pad,strel('disk',round(debrisrad),0));
    realnuc_la=bwlabel(DAs_pad);
    DAs_da=regionprops(realnuc_la,'Area','Centroid','PixelIdxList');
    %}
    [DAs_da,realnuc_la,finalcytoring]=buildcytoring(DAs_pad,REs_bs,nucr);
    fcr_da=regionprops(finalcytoring,'PixelIdxList');
    %%% screen out nuclei that are too small %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numcells=size(DAs_da,1);
    
    numrings=size(fcr_da,1);
    if numcells>numrings
        nucexclude=[numrings+1:numcells];
        DAs_da(nucexclude)=[];
    end
    nucexclude=zeros(numrings,1);
    for k=1:numrings
        if DAs_da(k).Area < minnucarea
            nucexclude(k)=1;
        end
    end
    nucexclude=find(nucexclude);
    DAs_da(nucexclude)=[];
    fcr_da(nucexclude)=[];
    numcells=size(DAs_da,1);
    %}
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    realnuc_la(ismember(realnuc_la,nucexclude))=0;
    finalcytoring(ismember(finalcytoring,nucexclude))=0;
    extractmask=bwmorph(realnuc_la,'remove') + logical(finalcytoring);
    %extractmask=bwmorph(realnuc_la,'remove');
    %extractmask=finalcytoring;
    %imwrite(uint16(extractmask),[cpdir,NEfile]);
    imshow(extractmask);
    keyboard;
    toc(timeframe)
end
toc(timetotal);