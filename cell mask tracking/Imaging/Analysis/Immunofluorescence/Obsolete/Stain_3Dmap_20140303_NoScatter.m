function Stain_3Dmap(conditions,datadir)
panel=0; %0:single image
if panel==0
    titlefont=8;
else
    titlefont=16;
end
dataset=20131216;
condnum=size(conditions,1);
%colorcode=distributecolormap(jet,condnum);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    Alldata=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                load([datadir,shot,'.mat'],'IFdata');
                Alldata=[Alldata;IFdata];
            end
        end
    end
    Hoechstval=Alldata(:,3).*Alldata(:,4);
    samplecells=Hoechstval>150000 & Hoechstval<250000;
    valuesx=log2(Alldata(samplecells,6));
    valuesy=log2(Alldata(samplecells,5));
    valuesz=Alldata(samplecells,7); valuesz(valuesz<1)=1; valuesz=log2(valuesz); valueszorg=valuesz;
    
    %hist(valuesz,1.95:0.05:6.95);xlim([2 7]);
    if dataset==20131216
        minz=0; anchorlowz=4; midz=7; anchorhighz=8.5; maxz=10; %serumrelease 20131126-pRb807811
    end
    
    polarizeoption=0;
    if polarizeoption==1
        removeid=find(valuesz>(midz-0.25) & valuesz<(midz+0.25));  %remove core
        %removeid=find(valuesz<(midz-0.25) | valuesz>(midz+0.25));  %remove poles
        valuesx(removeid)=[]; valuesy(removeid)=[]; valuesz(removeid)=[];
    end

    if panel==1
        switch dataset
            case 20131209
                subaxis(3,4,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02); %20131126
            case 20131216
                subaxis(3,3,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02); %20131209
            case 20140121  %20140121_p21CycA
                subaxis(3,4,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02); %20131209
        end
    end
    %%% select graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dscatter(valuesx,valuesy);
        %axis([17 19 5 11]); %p21
        %axis([17 19 2 12]); %CycA2
        %axis([17 19 3 10]); %CycD1
        %axis([17 19 5 8]); %CycE
    
    %scatter(valuesx,valuesy,20,valuesz,'fill'); colorbar;
    pretransit=valuesz> 1 & valuesz<5; posttransit=valuesz>5 & valuesz<8;
    selection=posttransit;
    scatter(valuesx(selection),valuesy(selection),40,valuesz(selection),'fill'); colorbar;
    
    %make3dmap(valuesx,valuesy,valuesz,minz,anchorlowz,midz,anchorhighz,maxz,valueszorg);
    %makedensitymap(valuesx,valuesy);
    %makeprobmap(valuesx,valuesy,valuesz,midz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if dataset==20131216
        hold on; line([0 15],[3 9],'linewidth',2,'color','r');  %20131126
        axis([2 10 2 7.5]); %20131126 Default:[4 9 3 9] NEDDi:[4.5 9.5 3.5 11.5];
        %axis([4.5 9.5 3 9]); %hr5 presentation
    end
    title(char(conditions(i,1)),'FontSize',titlefont);
end

if panel==0
    paperwidth=3; paperheight=3;
    xlabel('p21 log2(meanRFU)'); ylabel('CycD1 log2(meanRFU)');
else
    paperwidth=12; paperheight=10;
end
set(gcf,'color','w','PaperPosition',[0 0 paperwidth paperheight]); %6 5 (small panel: 4 3.3)
saveas(gcf,'h:\Downloads\Fig.jpg');
%{
figure,boxplot(axes,boxplotdata,offset,'labels',namecount);
title(titlestring);ylabel(descstring);
ylim([bmin bmax]);
set(gcf,'color','w','PaperPosition',[0 0 4 4]);
saveas(gcf,'h:\Downloads\Fig_boxplots.jpg');
%}
%{
figure;
errorbar(sigtime,sigmean,sigstd,'-mo','linewidth',2,'markeredgecolor','k','markerfacecolor',[.49 1 .63],'markersize',10);
title(titlestring);xlabel('time since serum release (hr)'); ylabel(descstring);
set(gcf,'color','w','PaperPosition',[0 0 4 4]);
saveas(gcf,'h:\Downloads\Fig_boxplots.jpg');
%}
end

function make3dmap(valuesx,valuesy,valuesz,minz,anchorlowz,midz,anchorhighz,maxz,orgvals)
rangez=maxz-minz;
normvaluesz=(valuesz-minz)/rangez;
cmapz=normvaluesz*64;
cmapz=round(cmapz); cmapz(cmapz<1)=1; cmapz(cmapz>64)=64;
lowcode=[0 0 1]; highcode=[1 1 0]; midcode=mean([lowcode;highcode]);
anchorlow=round(64*(anchorlowz-minz)/rangez)+1;
anchormid=round(64*(midz-minz)/rangez);
anchorhigh=round(64*(anchorhighz-minz)/rangez);
cmap=makecmap(lowcode,midcode,highcode,anchorlow,anchormid,anchorhigh);
cmap=jet;
colorsz=cmap(cmapz,:);
scatter(valuesx,valuesy,20,colorsz,'fill'); %default demo:10 panel:20

%{
binz=minz:rangez/100:maxz;
figure;
orgvals(orgvals<minz)=minz; orgvals(orgvals>maxz)=maxz;
binfill=[binz fliplr(binz)];
c1=[0 0 1]; colortemp=c1+.7*(1-c1);
histdat=hist(orgvals,binz)';
fill(binfill,[histdat;zeros(length(binz),1)],colortemp,'edgecolor',colortemp,'FaceAlpha', 0.4);
hold on;

valuesz(valuesz<minz)=[]; valuesz(valuesz>maxz)=[];
hist(valuesz,binz); xlim([minz maxz]);
title('Goat pRb(807/811) log2(meanRFU)');
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\Fig_Hist.jpg');

figure; imagesc(1:64); colormap(cmap);
set(gca,'XTick',[]); set(gca,'YTick',[]);
set(gcf,'color','w','PaperPosition',[0 0 4 0.5]);
saveas(gcf,'h:\Downloads\Fig_Colorbar.jpg');
%}
end

function makedensitymap(valuesx,valuesy)
stepsize=0.2; %default 0.2
halfstep=0.5*stepsize;
binx=min(valuesx)+halfstep:stepsize:max(valuesx)-halfstep;
biny=min(valuesy)+halfstep:stepsize:max(valuesy)-halfstep;
numbinx=length(binx); numbiny=length(biny);
countmap=zeros(numbiny,numbinx);
vbinx=zeros(numel(countmap),1); vbiny=vbinx; vbinz=vbinx;
total=numel(valuesx);
cc=0;
for i=1:numbinx
    samplex=find(valuesx>=binx(i)-halfstep & valuesx<binx(i)+halfstep);
    for j=1:numbiny
        cc=cc+1;
        vbinx(cc)=binx(i);
        vbiny(cc)=biny(j);
        if isempty(samplex)
            continue;
        end
        samplebinidx=valuesy(samplex)>=biny(j)-halfstep & valuesy(samplex)<biny(j)+halfstep;
        samplebin=samplex(samplebinidx);
        samplecount=numel(samplebin);
        if samplecount<=1
            samplecount=0;
        end
        countmap(j,i)=100*samplecount/total;
        %countmap(j,i)=samplecount;
        vbinz(cc)=100*samplecount/total;
    end
end
% imagesc(binx,biny,countmap);
% set(gca,'YDir','normal');
% cmap=jet; cmap(1,:)=[1 1 1];

maxperc=5;
%maxperc=max(countmap(:))
if max(countmap(:))>maxperc
    fprintf('countmap clipped\n');
end
countmap=countmap/maxperc;
countmap=countmap*64; countmap=round(countmap); countmap(countmap<1)=1; countmap(countmap>64)=64;
imagesc(binx,biny,countmap);
set(gca,'YDir','normal');
cmap=jet; cmap(1,:)=[1 1 1];
colormap(cmap);
%colorbar;
end

function makeprobmap(valuesx,valuesy,valuesz,midz)
stepsize=0.1; %default 0.2
halfstep=0.5*stepsize;
binx=min(valuesx)+halfstep:stepsize:max(valuesx)-halfstep;
biny=min(valuesy)+halfstep:stepsize:max(valuesy)-halfstep;
numbinx=length(binx); numbiny=length(biny);
countmap=ones(numbiny,numbinx)*-1;
vbinx=zeros(numel(countmap),1); vbiny=vbinx; vbinz=vbinx;
cc=0;
for i=1:numbinx
    samplex=find(valuesx>=binx(i)-halfstep & valuesx<binx(i)+halfstep);
    for j=1:numbiny
        cc=cc+1;
        vbinx(cc)=binx(i);
        vbiny(cc)=biny(j);
        if isempty(samplex)
            continue;
        end
        samplebinidx=valuesy(samplex)>=biny(j)-halfstep & valuesy(samplex)<biny(j)+halfstep;
        samplebin=samplex(samplebinidx);
        samplecount=numel(samplebin);
        samplepos=sum(valuesz(samplebin)>midz);
        if samplecount<1 %default <3
            continue;
        end
        countmap(j,i)=samplepos/samplecount;
    end
end
countmap=countmap*64; countmap=round(countmap);
countmap(countmap<2 & countmap>=0)=2;
countmap(countmap<0)=1; %anywhere there are no datapoints set to white
countmap(countmap>64)=64;
imagesc(binx,biny,countmap);
set(gca,'YDir','normal');
cmap=jet;
%cmap=ones(64,3); midoffset=12; cmap(32-midoffset:32+midoffset,:)=zeros(2*midoffset+1,3);
cmap(1,:)=[1 1 1];
colormap(cmap);
%cmapz=mat2gray(vbinz)*64; cmapz=round(cmapz); cmapz(cmapz<1)=1; cmapz(cmapz>64)=64;
%cmap=jet;
%cmap(1,:)=[1 1 1]; %0:white
%colorsz=cmap(cmapz,:);
%scatter(vbinx,vbiny,80,colorsz,'fill','square');
end