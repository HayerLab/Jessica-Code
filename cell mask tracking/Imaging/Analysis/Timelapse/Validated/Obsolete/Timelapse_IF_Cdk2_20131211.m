function Timelapse_IF_Cdk2(conditions,datadir)
condnum=size(conditions,1);
minlength=10;
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces1=[];
    alltraces2=[];
    allIF1=[];
    allIF2=[];
    allstats=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [traces1,traces2,IF1,IF2,stats]=main(datadir,shot,minlength);
                alltraces1=[alltraces1;traces1];
                alltraces2=[alltraces2;traces2];
                allIF1=[allIF1;IF1];
                allIF2=[allIF2;IF2];
                allstats=[allstats;stats];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
[numtraces,numframes]=size(alltraces1);
sampleid=(1:numtraces)';
detectrisetime=1;
if detectrisetime
    %%% detect onset of CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relative=0; %returns Cdk2start time relative to trace start
    [~,Cdk2start,~,badtraces]=getCdk2features_mother(alltraces1,allstats,minlength,relative);
    Cdk2start(isnan(Cdk2start))=0;
    sampleid(badtraces>0)=[];
end
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sorttraces=1;
if sorttraces
    [Cdk2start,idx]=sort(Cdk2start(sampleid));
    ordered=sampleid(idx);
    alltraces1=alltraces1(ordered,:);
    alltraces2=alltraces2(ordered,:);
    allIF1=allIF1(ordered);
    allIF2=allIF2(ordered);
    allstats=allstats(ordered,:);
end
%%% measure minimum CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltraces1,1);
minval=ones(numtraces,1)*NaN;
for i=1:numtraces
    minval(i)=min(smooth(alltraces1(i,allstats(i,1)+4:allstats(i,1)+minlength-1))); %default 19
end
%figure,hist(minval,50);title(['minimum DHB value 1-',num2str(minlength/framesperhr),'hrs after mitosis']);
%%% remove IF outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=allIF1;
IFdata(IFdata<1)=1;
IFdata=log2(IFdata);
%hist(IFdata(:),100);
%[mu1,sig1]=normfit(IFdata);
mu1=nanmean(IFdata(:));
sig1=nanstd(IFdata(:));
outliers=IFdata>(mu1+3*sig1) | IFdata<(mu1-3*sig1); % | IFdata<1
IFdata(outliers)=[];
alltraces1(outliers,:)=[];
alltraces2(outliers,:)=[];
Cdk2start(outliers)=[];
allstats(outliers,:)=[];
minval(outliers)=[];
%%% align data to immunostain time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltraces1,1);
alignedtraces1=ones(numtraces,numframes)*NaN;
alignedtraces2=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    alignedtraces1(i,numframes-allstats(i,3)+1:numframes)=alltraces1(i,allstats(i,1):allstats(i,2));
    alignedtraces2(i,numframes-allstats(i,3)+1:numframes)=alltraces2(i,allstats(i,1):allstats(i,2));
end
maxlength=prctile(allstats(:,3),90);
alignedtraces1=alignedtraces1(:,numframes-maxlength+1:numframes);
alignedtraces2=alignedtraces2(:,numframes-maxlength+1:numframes);
%%% display graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracedata=alignedtraces1;
time=(numframes-maxlength-1:numframes)/framesperhr;
sampleid=1:numtraces;

Cdk2inc=minval>0.5;
time_Cdk2inc=allstats(Cdk2inc,3);
time_Cdk2low=allstats(~Cdk2inc,3);
IF_Cdk2inc=IFdata(Cdk2inc);
IF_Cdk2low=IFdata(~Cdk2inc);

Cdk2activated=Cdk2start~=0;
sampleid_Cdk2act=sampleid(Cdk2activated);
time_Cdk2act=Cdk2start(Cdk2activated);
reltime_Cdk2act=time_Cdk2act-allstats(Cdk2activated,1);
IF_Cdk2act=IFdata(Cdk2activated);
IF_Cdk2nonact=IFdata(~Cdk2activated);
timesince_Cdk2act=numframes-time_Cdk2act;
heatmapgraph=0;
if heatmapgraph
    tracedata(tracedata>2.5)=2.5;
    tracedata(tracedata<0.3)=0.3;
    tracedata(isnan(tracedata))=0.3; %make invalid data black
    imagesc(time,sampleid,tracedata);
    % bluecode=[0 0 1]; blackcode=[0 0 0]; yellowcode=[1 1 0];
    % cmap=makecmap(bluecode,blackcode,yellowcode);
    % colormap(cmap);
    cmap=colormap(jet);
    cmap(1,:)=[0 0 0];
    colormap(cmap);
    %%% overlay markers at timepoints of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    addmarks=1;
    if addmarks
        hold on;
        scatter(time_Cdk2act/framesperhr,sampleid_Cdk2act,8,'r','*');
    end
    %%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xlabel('Time relative to Mitosis (hrs)'); ylabel('Individual traces');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]);
    saveas(gcf,'h:\Downloads\Fig1.jpg');
    %%% make IF panel to accompany heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure,imagesc(1,sampleid,IFdata);
    set(gca,'XTick',[],'YTick',[]);
    set(gcf,'color','w','PaperPosition',[0 0 1 8]);
    saveas(gcf,'h:\Downloads\Fig2.jpg');
end
G0durationhist=0;
if G0durationhist
    figure,hist(reltime_Cdk2act/framesperhr,100);
    title('G0 durations for cells entering G1');
end
timesincerisegraph=1;
if timesincerisegraph
    figure, hold on;
    binstep=1; binmin=minlength/framesperhr; binmax=20;
    binscattershade(timesince_Cdk2act/framesperhr,IF_Cdk2act,binstep,binmin,binmax,'b');
    %scatterandfit(timesince_Cdk2act/framesperhr,IF_Cdk2act,'b','r',1,binmin);
    title('Cyclin D1 levels vs time since Cdk2 activation');
    xlabel('Time since Cdk2 activation (hrs)'); ylabel('Cyclin D1 (log2(medianRFU))');
    xlim([binmin 20]);
    set(gcf,'color','w','PaperPosition',[0 0 3 3]); saveas(gcf,'h:\Downloads\Fig3.jpg');
end
actvsnonactgraph=0;
if actvsnonactgraph
    figure, hold on;
    data=[IF_Cdk2nonact;IF_Cdk2act];
    offset=[ones(numel(IF_Cdk2nonact),1);2*ones(numel(IF_Cdk2act),1)];
    boxplot(data,offset,'labels',{'Cdk2 activated','Cdk2 not activated'});
    title('Cyclin D1 levels with or without Cdk2 activation');
    ylabel('Cyclin D1 (log2(medianRFU))');
    set(gcf,'color','w','PaperPosition',[0 0 3 3]); saveas(gcf,'h:\Downloads\Fig4.jpg');
end
incvslowgraph=0;
if incvslowgraph
    figure, hold on;
    binstep=1; binmin=minlength/framesperhr; binmax=20;
    binscattershade(time_Cdk2inc/framesperhr,IF_Cdk2inc,binstep,binmin,binmax,'b');
    binscattershade(time_Cdk2low/framesperhr,IF_Cdk2low,binstep,binmin,binmax,'r');
    %scatterandfit(time_Cdk2inc/framesperhr,IF_Cdk2inc,'b','b',binstep,binmin);
    %scatterandfit(time_Cdk2low/framesperhr,IF_Cdk2low,'r','r',binstep,binmin);
    xlim([minlength/framesperhr 20]); ylim([7 12]);
    xlabel('Time since Mitosis (hrs)'); ylabel('Cyclin D1 (log2(medianRFU))');
    title('Cyclin D1 levels for Cdk2-inc vs Cdk2-low cells');
    set(gcf,'color','w','PaperPosition',[0 0 3 3]); saveas(gcf,'h:\Downloads\Fig5.jpg');
    figure, hold on;
    data=[IF_Cdk2inc;IF_Cdk2low];
    offset=[ones(numel(IF_Cdk2inc),1);2*ones(numel(IF_Cdk2low),1)];
    boxplot(data,offset,'labels',{'Cdk2-inc','Cdk2-low'});
    title('Cyclin D1 levels for Cdk2-inc vs Cdk2-low cells');
    ylabel('Cyclin D1 (log2(medianRFU))'); ylim([7 12]);
    set(gcf,'color','w','PaperPosition',[0 0 3 3]); saveas(gcf,'h:\Downloads\Fig6.jpg');
end
end


function [traces1,traces2,IF1,IF2,stats]=main(datadir,shot,minlength)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
load([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
%%% gate traces based on signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%samplecells=find(~isnan(IFdata(:,1))); %take all cells (not just daughters)
samplecells=find(~isnan(IFdata(:,1)) & ~isnan(genealogy)); %take only daughters
signals=tracedata(samplecells,:,5);
stats=tracestats(samplecells,:);
%hist(signals(:),0:10:2010); xlim([0 2000]); title('median nuclear DHB levels');
minthresh=75; %200 for Sabrina and Kyuho
[samplecells,~]=remove_short_neg_noisy(samplecells,signals,stats,minlength,minthresh);
%%% gate on another signal (optional) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signals=tracedata(samplecells,:,5);
% stats=tracestats(samplecells,:);
% minthresh=200;
% [samplecells,samplestats]=remove_short_neg_noisy(samplecells,signals,stats,minlength,minthresh);
%%% calculate signal of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces1=tracedata(samplecells,:,6)./tracedata(samplecells,:,5);
traces2=tracedata(samplecells,:,6);
stats=tracestats(samplecells,:);
%%% calculate IF values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF1=IFdata(samplecells,7); %ring of w4
IF2=IFdata(samplecells,7); %nuc of w5
end
