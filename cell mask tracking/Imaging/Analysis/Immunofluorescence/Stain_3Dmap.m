function Stain_3Dmap(conditions,datadir)
panel=0; %0:single image
if panel==0
    titlefont=8;
else
    titlefont=16;
end
dataset='20131216';
%condnum=size(conditions,1);
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);
%colorcode=distributecolormap(jet,condnum);
for i=1:uniquecondnum
    Alldata=[];
    condrow=find(ismember(conditions(:,1),uniquenames{i}));
    for c=condrow'
        rowmat=cell2mat(conditions(c,2));
        colmat=cell2mat(conditions(c,3));
        sitemat=cell2mat(conditions(c,4));
        for row=rowmat
            for col=colmat
                for site=sitemat
                    %shot=wellnum2str(row,col,site);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    load([datadir,shot,'.mat'],'IFdata');
                    Alldata=[Alldata;IFdata];
                end
            end
        end
    end
    Hoechstval=Alldata(:,3).*Alldata(:,4);
    samplecells=Hoechstval>100000 & Hoechstval<300000; %DMSO4hr:[150k 250k] %20131217:[100k 300k]
    valuesx=Alldata(:,6);
    posvalx=valuesx>=1; valuesx(valuesx<1)=1; valuesx=log2(valuesx);
    valuesy=Alldata(:,5);
    posvaly=valuesy>=1; valuesy(valuesy<1)=1; valuesy=log2(valuesy);
    valuesz=Alldata(:,7);
    posvalz=valuesz>=1; valuesz(valuesz<1)=1; valuesz=log2(valuesz); valueszorg=valuesz;
    posval=posvalx & posvaly & posvalz;
    valuesx=valuesx(samplecells & posval);
    valuesy=valuesy(samplecells & posval);
    valuesz=valuesz(samplecells & posval); valueszorg=valueszorg(samplecells & posval);
    
    %hist(valuesz,1.95:0.05:6.95);xlim([2 7]);
    switch dataset
        case '20131217'
            minz=7; anchorlowz=8; midz=9.25; anchorhighz=11; maxz=12;
            rangex=[5 10]; rangey=[5 8.5];
        case '20131216'
            %minz=0; anchorlowz=4; midz=7; anchorhighz=8.5; maxz=10;
            %minz=2; anchorlowz=3; midz=6.75; anchorhighz=9.5; maxz=10;
            %rangex=[3 10]; rangey=[2 7.5];
            minz=6; anchorlowz=6; midz=9.3; anchorhighz=12; maxz=12;
            rangex=[6 10]; rangey=[5 8];
        case 'CDKhysteresis'
            minz=0; anchorlowz=2; midz=7; anchorhighz=10; maxz=12;
        case 'cycDp21drugpanel'
            minz=0; anchorlowz=2; midz=7; anchorhighz=10; maxz=12;
        case 'IL2mutantsYT-'
            minz=4.5; anchorlowz=5; midz=6.5; anchorhighz=8.5; maxz=9;
    end
    
    %%% gate extreme values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    badvals=valuesz<minz | valuesz>maxz;
    valuesx=valuesx(~badvals);
    valuesy=valuesy(~badvals);
    valuesz=valuesz(~badvals); valueszorg=valueszorg(~badvals);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    polarizeoption=0;
    if polarizeoption==1
        removeid=find(valuesz>(midz-0.25) & valuesz<(midz+0.25));  %remove core
        %removeid=find(valuesz<(midz-0.25) | valuesz>(midz+0.25));  %remove poles
        valuesx(removeid)=[]; valuesy(removeid)=[]; valuesz(removeid)=[];
    end

    if panel==1
        switch dataset
            case '20131209'
                subaxis(3,4,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02); %20131126
            case '20131216'
                %subaxis(3,3,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02);
                subaxis(2,4,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02,'SV',0.1);
            case '20140121'  %20140121_p21CycA
                subaxis(3,4,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02); %20131209
            case 'CDKhysteresis'
                subaxis(2,3,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02);
            case 'cycDp21drugpanel'
                subaxis(3,4,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02);
            case 'IL2mutantsYT-'
                subaxis(2,3,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02);
        end
    end
    %%% select graph %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dscatter(valuesx,valuesy);
        %axis([17 19 5 11]); %p21
        %axis([17 19 2 12]); %CycA2
        %axis([17 19 3 10]); %CycD1
        %axis([17 19 5 8]); %CycE
    
    %scatter(valuesx,valuesy,20,valuesz,'fill');
    %pretransit=valuesz> 1 & valuesz<5; posttransit=valuesz>5 & valuesz<8;
    %selection=posttransit;
    %scatter(valuesx(selection),valuesy(selection),40,valuesz(selection),'fill'); colorbar;
    
    
    %hist(valuesz,50); hold on; line([midz midz],[0 54],'color','r','linewidth',4); xlabel('pRb log2(meanRFU)'); ylabel('# cells');
    stoich=valuesy./valuesx;
    
    binarytransform(stoich,valuesz,midz);
    %make3dmap(valuesx,valuesy,valuesz,minz,anchorlowz,midz,anchorhighz,maxz,valueszorg,rangex,rangey);
    %makedensitymap(valuesx,valuesy);
    %makeprobmap(valuesx,valuesy,valuesz,midz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hold on;
    switch dataset
        case '20131216'
            line([0 15],[3 9],'linewidth',2,'color','r');
            axis([2 10 2 7.5]); %20131126 Default:[4 9 3 9] NEDDi:[4.5 9.5 3.5 11.5];
        case 'CDKhysteresis'
            line([0 10],[1.6 3.1],'linewidth',2,'color','r');
            axis([0 9 1.6 3.0]);
        case 'cycDp21drugpanel'
            line([0 12],[2.5 9],'linewidth',2,'color','r');
            axis([0 9 2 9]);
        case 'IL2mutantsYT-'
            %line([0 12],[2.5 9],'linewidth',2,'color','r');
            axis([6 12 4 9]);
    end
    title(char(uniquenames{i}),'FontSize',titlefont);
end

if panel==0
    paperwidth=3; paperheight=3;
    %xlabel('p21 log2(meanRFU)'); ylabel('CycD1 log2(meanRFU)');
    xlabel('pERK log2(meanRFU)'); ylabel('pAKT log2(meanRFU)');
else
    paperwidth=12; paperheight=8; %[12 10]
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

function binarytransform(input,output,midpoint)
numbins=50;
bmin=min(input); bmax=max(input); bstep=(bmax-bmin)/numbins;
bin=bmin:bstep:bmax;
bcalc=ones(numbins,1)*NaN;
for i=1:numbins
    inputidx=input>bin(i) & input<bin(i+1);
    ovals=output(inputidx);
    bcalc(i)=sum(ovals>midpoint)/length(ovals);
end
truebin=bin+bstep/2;
scatter(truebin(1:numbins),bcalc,100,'markeredgecolor','k','markerfacecolor',[.49 1 .63]);
xlim([0.5 1.2]);
xlabel('cycD1/p21 stoichiometry');
ylabel('% pRb-positive')
set(gcf,'color','w','PaperPosition',[0 0 4 4]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end

function make3dmap(valuesx,valuesy,valuesz,minz,anchorlowz,midz,anchorhighz,maxz,orgvals,rangex,rangey)
rangez=maxz-minz;
normvaluesz=(valuesz-minz)/rangez;
cmapz=normvaluesz*64;
cmapz=round(cmapz); cmapz(cmapz<1)=1; cmapz(cmapz>64)=64;
lowcode=[0 0 1]; highcode=[1 1 0]; midcode=mean([lowcode;highcode]);
anchorlow=round(64*(anchorlowz-minz)/rangez)+1;
anchormid=round(64*(midz-minz)/rangez);
anchorhigh=round(64*(anchorhighz-minz)/rangez);
cmap=makecmap(lowcode,midcode,highcode,anchorlow,anchormid,anchorhigh);
%cmap=jet;
colorsz=cmap(cmapz,:);
scatter(valuesx,valuesy,20,colorsz,'fill'); %default demo:10 panel:20
axis([rangex rangey]);
keyboard;
%{
%%% plot p21 vs CycD
scatter(valuesx,valuesy,20,'o','fill');
%axis([3 10 3.5 7.5]);
xlabel('p21 log2(meanRFU)');
ylabel('CycD1 log2(meanRFU)');
set(gcf,'color','w','PaperPosition',[0 0 3 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');

%%% histogram
binz=minz:rangez/40:maxz;
hist(valuesz,binz); xlim([minz maxz]);
title('pRb(807/811) log2(meanRFU)');
set(gcf,'color','w','PaperPosition',[0 0 4 3]);
saveas(gcf,'h:\Downloads\Fig_Hist.jpg');
%%% colorbar
figure; imagesc(1:64); colormap(cmap);
set(gca,'XTick',[]); set(gca,'YTick',[]);
set(gcf,'color','w','PaperPosition',[0 0 4 0.5]);
saveas(gcf,'h:\Downloads\Fig_Colorbar.jpg');

binz=minz:rangez/40:maxz;
figure;
orgvals(orgvals<minz)=minz; orgvals(orgvals>maxz)=maxz;
binfill=[binz fliplr(binz)];
c1=[0 0 1]; colortemp=c1+.7*(1-c1);
histdat=hist(orgvals,binz)';
fill(binfill,[histdat;zeros(length(binz),1)],colortemp,'edgecolor',colortemp,'FaceAlpha', 0.4);
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