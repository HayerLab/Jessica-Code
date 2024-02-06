function sample_measure()
cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timelapse\hDHB Optimization\hDHB_Geminin_20x_Steve\';  %folder containing the movies folders [CHANGE]
datadir = ([path,'Data\']);
cpdir = [path,'CroppedProcessed\'];
nucr=16; %MCF10A/10x:8 MCF10A/20x:16 HS68andHeLa/10x:9
movie=1;
singlecell=0;
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
%for f=framezero+106
%for f=framezero+108
%for f=framezero+117
for f=framezero+77
    %% read in pre-cropped images
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
    REs_bs=bgsub(REs_or,10*nucr,0.05);
    %CEs_bl=(imfilter(CEs_or,fspecial('disk',floor(nucr/2)),'symmetric'));
    %CEs_bs=bgsub(CEs_bl,10*nucr,0.05);
    
    %% get data
    %DAs_ma=getdapimask(DAs_bs,nucr);  %make the mask. this segmentation function is the engine of Feng Chiao's script
    DAs_ma=getnucmask(DAs_bs,nucr); %Mingyu's segmentation script
    
    %{
    %%% watershed mask
    D=bwdist(~DAs_ma);
    D=-D;
    D(~DAs_ma)=-Inf;
    DAs_ws=watershed(D);
    DAs_ma=DAs_ws & DAs_ma;
    %}
    
    %DAs_la=bwlabel(DAs_ma);  %labels the objects with numbers
    
    if singlecell
        centralnuc=DAs_la==DAs_la(floor(height/2),floor(width/2));
        DAs_la=bwlabel(centralnuc);
    end
    
    %orgring=imdilate(DAs_la,strel('disk',4,8)) - imdilate(DAs_la,strel('disk',1,8));
    %BAs_la = imdilate(DAs_la,strel('disk',2,8));
    %MAs_la=imerode(DAs_la,strel('disk',4,8));
    %newring=DAs_la-imerode(DAs_la,strel('disk',3,8)); %corrected to match width of orgring
    %newring=BAs_la-imerode(BAs_la,strel('disk',3,8)); %corrected to match width of orgring
    
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
        offset=round(nucr/8);
        orderedsetoffset=[orderedset(offset+1:end,:);orderedset(1:offset,:)];
        diffset=orderedsetoffset-orderedset;
        gradset=atan2(diffset(:,2),diffset(:,1));   %angle in radians
        gradsetoffset=[gradset(offset+1:end,:);gradset(1:offset,:)];
        graddiffset=gradsetoffset-gradset;
        graddiffset=graddiffset+2*pi*(graddiffset<0);
        gradthresh = pi/3;
        %inversegradthresh = -pi+gradthresh-pi;
        %vIdx=(graddiffset>gradthresh & graddiffset<pi) | (graddiffset>inversegradthresh & graddiffset<-pi);
        vIdx=graddiffset>gradthresh & graddiffset<pi;
        %{
        vposvec=orderedset(vIdx+offset,:);
        for vc=1:size(vposvec,1)
            verticesmask(vposvec(vc,2),vposvec(vc,1))=1;
        end
        %}
        %%% find center of contiguous points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %consider vIdxobs=imclose(vIdxobs);
        vIdxobs=regionprops(logical(vIdx),'Centroid');
        vIdx=zeros(length(vIdxobs),1);
        for rpc=1:length(vIdxobs)
            vIdx(rpc)=floor(vIdxobs(rpc).Centroid(2));
        end
        %%% if odd # vertices, remove lowest angle change %%%%%%%%%%%%%%%%%
        if rem(length(vIdx),2)                      %odd number of vertices
            vgraddiff=graddiffset(vIdx);
            %{
            for vgi=1:length(vgraddiff)
                if vgraddiff(vgi)<0
                    vgraddiff(vgi)=vgraddiff(vgi)+2*pi;
                end
            end
            %}
            %[~,ordx]=sort(vgraddiff);
            [~,idx]=min(vgraddiff);
            %vIdx=vIdx(ordx);
            %vIdx(1)=[];
            vIdx(idx)=[];
        end
        %%% pair and connect vertices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if length(vIdx)>=2
            while length(vIdx)>=2
                vpos=orderedset(vIdx+offset,:);
                vposadj=[vpos(2:end,:);vpos(1,:)];
                vgrad=gradset(vIdx);
                pair=vposadj-vpos;
                pairdist=sum(pair.^2,2);
                pairgrad=atan2(pair(:,2),pair(:,1));
                pairgraddiff=pairgrad-vgrad;
                
                pairgraddiff=pairgraddiff+2*pi*(pairgraddiff<0);
                vlo=pairgraddiff>gradthresh & pairgraddiff<pi;
                
                [~,ordx]=sort(pairdist);
                vlo=vlo(ordx);
                idx=min(ordx(~vlo));
                %[~,idx]=min(pairdist);
                idxadj=idx+1;
                if idxadj>length(vIdx)
                    idxadj=1;
                end
                [bcx,bcy]=bridge(vpos(idx,:),vpos(idxadj,:));
                for bci=1:length(bcx)
                    verticesmask(bcy(bci),bcx(bci))=1;
                end
                vIdx([idx,idxadj])=[];
            end
            %{
            vpos=orderedset(vIdx+offset,:);
            [bcx,bcy]=bridge(vpos(1,:),vpos(2,:));
            for bci=1:length(bcx)
                verticesmask(bcy(bci),bcx(bci))=1;
            end
            %}
            
            %{
            lengthx=vpos(2,1)-vpos(1,1);
            lengthy=vpos(2,2)-vpos(1,2);
            longerside=max([abs(lengthx) abs(lengthy)]);
            stepx=lengthx/longerside;
            stepy=lengthy/longerside;
            for bs=0:longerside
                bx=vpos(1,1)+round(bs*stepx);
                by=vpos(1,2)+round(bs*stepy);
                verticesmask(by,bx)=1;
            end
            %}
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    orgmoviering=(orgmoviering>0)*10;
    orgmoviering=orgmoviering(2:end-1,2:end-1);
    verticesmask=bwmorph(verticesmask,'dilate',1);
    verticesmask=(verticesmask>0)*10;
    verticesmask=verticesmask(2:end-1,2:end-1);
    %verticesmask=bwlabel(verticesmask);
    %verticesmask=imdilate(verticesmask,strel(4,8));
    
    
    %% Add frame to movie
    if movie
        
        %tempframe=imadjust(mat2gray(DAs_bs));
        tempframe=imadjust(mat2gray(DAs_bs));
        tempframe(:,:,2)=imadjust(mat2gray(orgmoviering));
        tempframe(:,:,3)=imadjust(mat2gray(verticesmask));
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
end
%%
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