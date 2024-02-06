function sample_measure()
cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_20x_Steve\';  %folder containing the movies folders [CHANGE]
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
nucr=16; %MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
movie=1;
bridgeimage=0;
%%%% sample specification %%%%%%%%%%%
row=2;
col=1;
site=1;
track=779;
%%% setup tempwell %%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
shottrack=[shot,'_track',num2str(track)];
load([datadir,shot,'_alldata'],'bestsp','best_rc','corruptlist','leaveoutlist');
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
framezero = 63;

if movie
    M=VideoWriter([datadir,shottrack],'Uncompressed AVI');
    M.FrameRate = 4;
    open(M);
end
tempread=imread([cpdir,shot,'_nucedge_',num2str(1),'.tif']);
[m,n]=size(tempread);

%for f=best_rc(track,1):best_rc(track,3)
%for f=framezero+143
%for f=framezero+90
for f=framezero+152
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
    [height,width]=size(DAs_or);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% image processing
    %DAs_bl=log(imfilter(DAs_or,fspecial('disk',floor(nucr/2)),'symmetric')); 
    DAs_bs=bgsub_test(DAs_or,10*nucr,0.05);  
    %REs_bl=(imfilter(REs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    REs_bs=bgsub_test(REs_or,10*nucr,0.05);
    %CEs_bl=(imfilter(CEs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    %CEs_bs=bgsub(CEs_bl,10*nucr,0.05);
    
    %% get data
    %DAs_ma=getdapimask(DAs_bs,nucr);  %make the mask. this segmentation function is the engine of Feng Chiao's script
    DAs_ma=getnucmask(DAs_bs,nucr); %Mingyu's segmentation script
    REs_ma=ThreshImage_test(REs_bs,0);  %looser threshold
    %REs_ma=single(REs_bs>0);
    REs_ma=imopen(REs_ma,strel('disk',nucr/4,8));
    REs_ma=bwmorph(REs_ma,'remove');

    %orgmoviering=bwmorph(DAs_la,'remove');
    DAs_pad=zeros(size(DAs_ma,1)+2,size(DAs_ma,2)+2);
    DAs_pad(2:end-1,2:end-1)=DAs_ma;
    DAs_pad=imopen(DAs_pad,strel('disk',nucr/2));  %this would remove stuff smaller than 25% of normal cell size.  Also removes spurs at edges & protruding points (but not intruding points, I think).
    orgmoviering=bwmorph(DAs_pad,'remove');
    [ringlabels,obnum]=bwlabel(orgmoviering);
    verticesmask=zeros(size(orgmoviering));
    
    %%% vertice detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ci=1:obnum
        [r,c]=find(ringlabels==ci);
        set=[c,r]; %adjust to x-y convention
        orderedset=zeros(size(set));
        orderedset(1,:) = set(1,:);
        numpoints=length(r);
        %%% determine clockwise direction from r(1),c(1) %%%%%%%%%%%%%%%%%%
        x=set(1,1); y=set(1,2);
        %adj=[x+1,y;x+1,y+1;x,y+1;x-1,y+1;x-1,y;x-1,y-1;x,y-1;x+1,y-1]; %clockwise from pos x axis
        adj=[x+1,y;x+1,y-1;x,y-1;x-1,y-1;x-1,y;x-1,y+1;x,y+1;x+1,y+1]; %clockwise from pos x axis (wrt y-axis convention)
        [Pradj,Ixset]=ismember(adj,set,'rows'); %present in adj, index in set
        if sum(Pradj)~=2
            fprintf(['Other than 2 neighbors detected!\n']);
        end
        poles=find(Pradj);
        side=round(mean(poles));
        nucpres = DAs_pad(adj(side,2),adj(side,1));
        if nucpres
            pole = poles(1);
        else
            pole = poles(2);
        end
        current=set(Ixset(pole),:);  %update
        set(1,:)=[];                 %remove previous center
        %%% order ring coordinates contiguously %%%%%%%%%%%%%%%%%%%%%%%%%%%
        for pp=2:numpoints-1
            set((set(:,1)==current(1) & set(:,2)==current(2)),:)=[];    %remove previous center
            orderedset(pp,:)=current;
            x=current(1); y=current(2);
            adj=[x+1,y;x+1,y+1;x,y+1;x-1,y+1;x-1,y;x-1,y-1;x,y-1;x+1,y-1]; %clockwise from pos x axis
            [Pradj,Ixset]=ismember(adj,set,'rows'); %present in adj, index in set
            if sum(Pradj)~=1
                fprintf(['Other than 1 neighbor detected!\n']);
            end
            pole=find(Pradj);
            current=set(Ixset(pole),:);  %update
        end
        orderedset(numpoints,:)=current;
        %%% calculate tangent angle & detect vertices %%%%%%%%%%%%%%%%%%%%%
        offsetshort=round(nucr/4);
        offsetlong=2*offsetshort;
        halfoffset=0; %was round(offset/2)
        gradientoffsetshort=offsetshort; %was offset*2
        gradientoffsetlong=offsetlong;
        
        orderedsetoffsetshort=[orderedset(offsetshort+1:end,:);orderedset(1:offsetshort,:)];
        orderedsetoffsetlong=[orderedset(offsetlong+1:end,:);orderedset(1:offsetlong,:)];
        shortdiff=orderedsetoffsetshort-orderedset;
        longdiff=orderedsetoffsetlong-orderedset;
        shortgrad=atan2(shortdiff(:,2),shortdiff(:,1));   %angle in radians
        longgrad=atan2(longdiff(:,2),longdiff(:,1));
        shortgradoffset=[shortgrad(gradientoffsetshort+1:end,:);shortgrad(1:gradientoffsetshort,:)];
        longgradoffset=[longgrad(gradientoffsetlong+1:end,:);longgrad(1:gradientoffsetlong,:)];
        shortgraddiff=shortgradoffset-shortgrad;
        longgraddiff=longgradoffset-longgrad;
        shortgraddiff=shortgraddiff+2*pi*(shortgraddiff<0);
        longgraddiff=longgraddiff+2*pi*(longgraddiff<0);
        shortgradthresh = pi/6; %was pi/4
        longgradthresh = pi/6;
        vIdxmasklong=longgraddiff>longgradthresh & longgraddiff<pi;
        vIdxmasklong=[zeros(offsetlong,1);vIdxmasklong(1:end-offsetlong)];
        vIdxmasklong=imdilate(vIdxmasklong,strel('square',1+nucr));
        shortgraddiff(shortgraddiff>=pi)=0;
        %vIdxmaskshort=shortgraddiff>shortgradthresh & shortgraddiff<pi;
        vIdxmaskshort=shortgraddiff>shortgradthresh;
        %vIdxmaskshort=getmax(vIdxmaskshort,shortgraddiff);
        
        vIdxmaskshort=imclose(vIdxmaskshort,strel('square',3));
        vIdxobs=regionprops(bwlabel(vIdxmaskshort),'PixelIdxList');
        maxmask=zeros(size(vIdxmaskshort));
        %vIdx=zeros(length(vIdxobs),1);
        for rpc=1:length(vIdxobs)
            %vIdx(rpc)=floor(vIdxobs(rpc).Centroid(2));
            pix=vIdxobs(rpc).PixelIdxList;
            [~,index]=max(shortgraddiff(pix));
            maxmask(pix(index)+offsetshort)=1;
        end
        
        vIdxmask=vIdxmasklong & maxmask;
        vIdx=find(vIdxmask);
        %%% mark candidate vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if bridgeimage
            vpos=orderedset(vIdx,:);
            for vc=1:size(vpos,1)
                verticesmask(vpos(vc,2),vpos(vc,1))=1;
            end
        end
        %%% pair and connect vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if length(vIdx)>=2
            periIdx=vIdx;
            perisize=length(vIdxmask);
            periIdxadj=[periIdx(2:end);perisize+periIdx(1)];
            pairperi=periIdxadj-periIdx;
        end
        while length(vIdx)>=2
            vpos=orderedset(vIdx,:);
            vposadj=[vpos(2:end,:);vpos(1,:)];
            pair=vposadj-vpos;
            pairdist=sqrt(sum(pair.^2,2));
            
            [~,ordx]=sort(pairperi./pairdist);
            idx=ordx(end);
            
            idxadj=idx+1;
            if idxadj>length(vIdx)
                idxadj=1;
            end
            
            if length(vIdx)<=3
                if length(vIdx)==3
                    idxadjadj=idxadj+1;
                    if idxadjadj>length(vIdx)
                        idxadjadj=1;
                    end
                    pairperi(idxadj)=pairperi(idxadj)+pairperi(idxadjadj);
                    pairperi(idxadjadj)=nucr*pi;
                end
                vlo=pairperi<nucr*pi; %less than half circumference of avg size nucleus
                if sum(vlo)         %either side too short
                    break
                end
            end
            
            [bcx,bcy]=bridge(vpos(idx,:),vpos(idxadj,:));
            for bci=1:length(bcx)
                verticesmask(bcy(bci),bcx(bci))=1;
            end
            
            previdx=idx-1;
            if previdx==0
                previdx=length(vIdx);
            end
            pairperi(previdx)=pairperi(previdx)+length(bcx)-1+pairperi(idxadj);
            
            vIdx([idx,idxadj])=[];
            pairperi([idx,idxadj])=[];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DAs_pad=DAs_pad(2:end-1,2:end-1);
    verticesmask=imdilate(verticesmask,strel('disk',1,8)); %required because bwlabel can't discriminate a one-pixel separation
    verticesmask=verticesmask(2:end-1,2:end-1);
    if bridgeimage
        orgmoviering=(orgmoviering>0)*10;
        orgmoviering=orgmoviering(2:end-1,2:end-1);
    end
    if ~bridgeimage
        DAs_pad=DAs_pad & ~verticesmask;
        orgmoviering=bwmorph(DAs_pad,'remove');
    end
    
    outerrad=nucr/4;
    cytoringoutermass=bwmorph(DAs_pad,'thicken',outerrad);
    cytoringoutermass(1,:)=0;cytoringoutermass(end,:)=0;cytoringoutermass(:,1)=0;cytoringoutermass(:,end)=0;
    cro_la=bwlabel(cytoringoutermass);
    cri_la=imerode(cro_la,strel('disk',outerrad/2+1,0));     %change this back before running!
    cytonuc_la=imerode(cro_la,strel('disk',1.5*outerrad,0)); %make sure calc inside true nuc
    realnuc_la=cro_la.*DAs_pad;
    cytoring=cro_la-cri_la;
    
    borders=bwmorph(cytoringoutermass,'bothat');
    borders=imdilate(borders,strel('disk',outerrad,8));
    cytoring=cytoring.*(~borders);
    
    %%% only for visualization %%%%%%%%%%%%%%%%%
    cytoringouter=bwmorph(cro_la,'remove');
    %cytoringouter=cytoringouter.*(~borders);
    %cytoringinner=bwmorph(cri_la,'remove');
    %cytoringinner=cytoringinner.*(~borders);
    cytonuc=bwmorph(cytonuc_la,'remove');
    realnuc=bwmorph(realnuc_la,'remove');
    %cytoringedges=cytoringinner | cytoringouter | cytonuc;
    %cytoringedges=cytoringinner | cytoringouter;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DAs_da=regionprops(cytonuc_la,'Area','Centroid','PixelIdxList'); %finds the centroid,etc of each labeled object  %to test, type "DAs_da(1).Area"
    ringxypos=regionprops(cytoring, 'PixelIdxList','PixelList');   %PixelList=[coords from left, coords from top]
    
    %%% find spoke edges %%%%%%%%
    numcells=size(DAs_da,1);
    inspokelength=nucr/2;
    outspokelength=nucr*2;
    tmap=zeros(height,width);%tbmap=zeros(height+outspokelength,width+outspokelength);
    E_home=tmap;E_in=tmap;%E_out=tbmap;
    for i=1:numcells
        cenx=round(DAs_da(i).Centroid(1));  %Centroid=[coords from left edge, coords from top]
        ceny=round(DAs_da(i).Centroid(2));
        ringcoords=ringxypos(i).PixelList;
        EWcandidates=ringcoords(ringcoords(:,2)==ceny,1);
        ty=ceny;tx=max(EWcandidates);
        E_home(ty,tx)=i; E_in(ty,tx-inspokelength:tx)=i; %E_out(ty,tx:tx+outspokelength)=i;
    end
    nhood=zeros(1,2*outspokelength+1);nhood(outspokelength+2:end)=1;
    E_out=imdilate(E_home,strel(nhood));
    %E_out=E_out(1:height,1:width);
        
    cytogradients=edge(REs_bs,'canny');
    cytogradients=bwmorph(cytogradients,'thin',Inf);
    %cytogradients=bwmorph(cytogradients,'diag'); %removes 8-connectivity of background
    cg_da=regionprops(cytogradients,'BoundingBox','PixelIdxList');
    for i=1:size(cg_da,1)
        longestbound=max(cg_da(i).BoundingBox([3 4]));
        if longestbound < nucr*2
            cytogradients(cg_da(i).PixelIdxList)=0;
        end
    end
    
    %%% remove ring outside of cytogradients %%%%%%%%%
    boundthresh=nucr/4;
    legitedges=bwlabel(cytogradients);
    %numedges=max(max(legitedges));
    
    %{
    innerscreen=legitedges.*logical(cytonuc_la);
    nucgradients=unique(innerscreen);
    for j=1:length(nucgradients)
        legitedges(legitedges==nucgradients(j))=0;
    end
    %}
    
    %{
    innerscreen=legitedges.*cytoringoutermass;
    pseudonuc=regionprops(innerscreen,'Area');
    highestlabel=max(max(innerscreen));
    for j=1:highestlabel
        if pseudonuc(j).Area > nucr*2*pi
            legitedges(legitedges==j)=0;
        end
    end
    %}
    
    gradientratio=4;
    innerscreen=legitedges.*cytoringoutermass;
    outerscreen=legitedges.*~cytoringoutermass;
    pseudonuc=regionprops(innerscreen,'Area');
    pseudocyto=regionprops(outerscreen,'Area');
    highestlabel=max(max(innerscreen));
    for j=1:highestlabel
        if pseudonuc(j).Area > gradientratio*pseudocyto(j).Area
            legitedges(legitedges==j)=0;
        end
    end
    
    inbounds=cytoringoutermass & ~legitedges & ~REs_ma;
    inbounds=imopen(inbounds,strel('disk',boundthresh,8));
    
    %cytoringouter=cytoringouter.*inbounds;
    cytoring=cytoring.*inbounds;
    
    %%% only for visualization %%%%%%%%%%%%%%%%%
    %cytoringperi=bwmorph(cytoring,'remove');
    cytoringperi=cytoring;
    %cytoringperi=cytoringouter;
    cytoringedges=logical(legitedges) | REs_ma;
    %cytoringedges=cytogradients;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{
    mz = size(DAs_da,1);
    XX=zeros(mz,1);YY=zeros(mz,1);  %define vectors which will hold the info from the structured array, DAs_da 
    AC=zeros(mz,1);%PP=zeros(size(DAs_da,1),1);DI=zeros(size(DAs_da,1),1);
    DD=zeros(mz,1);RR=zeros(mz,1); %FR=zeros(size(DAs_da,1),1); %DD is a vector of median dapi intensities for each object in log scale
    avgringyfp=zeros(mz,1);
    for cc=1:mz   %run a loop to put the info from the structured array in to the vectors; this is for each cell
        XX(cc,1)=DAs_da(cc).Centroid(1);  %x value of centroid
        YY(cc,1)=DAs_da(cc).Centroid(2);  %y value of centroid
        AC(cc,1)=DAs_da(cc).Area;
        DD(cc,1)=median(DAs_bs(DAs_da(cc).PixelIdxList));
        RR(cc,1)=median(REs_bs(DAs_da(cc).PixelIdxList));  %nuc yfp intensity

        allringpixels=REs_bs(ringxypos(cc).PixelIdxList);
        topringpixels=allringpixels(allringpixels>=prctile(allringpixels, 50));  %get top 50th percentile of ring pixels                
        avgringyfp(cc,1)=mean(topringpixels);  %take mean of pixels in the top 50th pctile
    end
    %}
    
    %% Add frame to movie
    if movie
        
        %tempframe=imadjust(mat2gray(DAs_bs));
        %tempframe=imadjust(mat2gray(DAs_bs));
        %tempframe(:,:,2)=imadjust(mat2gray(orgmoviering));
        %tempframe(:,:,3)=imadjust(mat2gray(verticesmask));
        %tempframe(:,:,3)=imadjust(mat2gray(REs_bs));
        %tempframe(:,:,3)=0;
        %}
        
        
        %DAs_bs=imfilter(DAs_bs,fspecial('disk',1),'symmetric');
        heatimage = imadjust(mat2gray(REs_bs));
        normimage = heatimage/(max(max(heatimage)));
        codeimage = round(normimage*128);
        codeimage(codeimage==0) = 1;
        
        blackring=1;
        %blackimage=cytoringedges | cytogradients;
        blackimage=cytoringedges | realnuc;
        whiteimage=cytoringperi;
        if ~blackring
            tempframe = colorcode1(codeimage);
            tempframe(:,:,2) = colorcode2(codeimage);
            tempframe(:,:,3) = colorcode3(codeimage);
        else
            tempframe = colorcode1(codeimage);
            tempframe(logical(blackimage))=0;
            tempframe(logical(whiteimage))=1;
            frame2=colorcode2(codeimage);
            frame2(logical(blackimage))=0;
            frame2(logical(whiteimage))=0;
            tempframe(:,:,2) = frame2;
            frame3=colorcode3(codeimage);
            frame3(logical(blackimage))=0;
            frame3(logical(whiteimage))=0;
            tempframe(:,:,3) = frame3;
        end
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
    %set(gcf,'PaperPosition',[0 0 3 3]);
    %saveas(gcf,'h:\Downloads\Fig.jpg');
end
cd ..\Analysis; %return to this directory

end

function [bx,by] = bridge(vpos1,vpos2)
lengthx=vpos2(1)-vpos1(1);
lengthy=vpos2(2)-vpos1(2);
longerside=max([abs(lengthx) abs(lengthy)]);
stepx=lengthx/longerside;
stepy=lengthy/longerside;
bx=zeros(longerside+1,1);by=zeros(longerside+1,1);
for bs=0:longerside
    bx(1+bs)=vpos1(1)+round(bs*stepx);
    by(1+bs)=vpos1(2)+round(bs*stepy);
end
end