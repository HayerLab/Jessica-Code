function Timelapse_IF_Cdk2(conditions,datadir)
condnum=size(conditions,1);
minlength=15;
IFstring='CycD1';
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
                shot=wellnum2str(row,col,site);
                %shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
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
    [dummy,Cdk2start,dummy,badtraces]=getCdk2features(alltraces1,allstats,minlength,relative);
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
%%% remove IF outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=allIF1;
%hist(IFdata(:),100); %check this before performing log2
IFdata(IFdata<1)=1;
IFdata=log2(IFdata);
%hist(IFdata(:),100);
mu1=nanmean(IFdata(:));
iqr1=prctile(IFdata(:),75)-prctile(IFdata(:),25);
outliers=IFdata>mu1+3*iqr1 | IFdata<mu1-3*iqr1;

IFdata(outliers)=[];
alltraces1(outliers,:)=[];
alltraces2(outliers,:)=[];
Cdk2start(outliers)=[];
allstats(outliers,:)=[];
%%% align data to immunostain time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltraces1,1);
IFalignedtraces1=ones(numtraces,numframes)*NaN;
IFalignedtraces2=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    IFalignedtraces1(i,numframes-allstats(i,3)+1:numframes)=alltraces1(i,allstats(i,1):allstats(i,2));
    IFalignedtraces2(i,numframes-allstats(i,3)+1:numframes)=alltraces2(i,allstats(i,1):allstats(i,2));
end
maxlength=prctile(allstats(:,3),90);
IFalignedtraces1=IFalignedtraces1(:,numframes-maxlength+1:numframes);
IFalignedtraces2=IFalignedtraces2(:,numframes-maxlength+1:numframes);
%%% display graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces_IFaligned=IFalignedtraces1;
time_IFaligned=(numframes-maxlength-1:numframes)/framesperhr;
sampleid=1:numtraces;

Cdk2activated=Cdk2start~=0;
sampleid_Cdk2act=sampleid(Cdk2activated);
time_Cdk2act=Cdk2start(Cdk2activated);
allstats_Cdk2act=allstats(Cdk2activated,:);
IF_Cdk2act=IFdata(Cdk2activated);
IF_Cdk2nonact=IFdata(~Cdk2activated);
timesince_Cdk2act=numframes-time_Cdk2act;
heatmapgraph=0;
if heatmapgraph
    traces_IFaligned(traces_IFaligned>2.5)=2.5;
    traces_IFaligned(traces_IFaligned<0.3)=0.3;
    traces_IFaligned(isnan(traces_IFaligned))=0.3; %make invalid data black
    imagesc(time_IFaligned,sampleid,traces_IFaligned);
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
    %saveas(gcf,'h:\Downloads\Fig1.jpg');
    %%% make IF panel to accompany heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure,imagesc(1,sampleid,IFdata);
    set(gca,'XTick',[],'YTick',[]);
    set(gcf,'color','w','PaperPosition',[0 0 1 8]);
    %saveas(gcf,'h:\Downloads\Fig2.jpg');
end
daughters_Cdk2act=~isnan(allstats(Cdk2activated,4));
reltime_Cdk2act=time_Cdk2act(daughters_Cdk2act)-allstats_Cdk2act(daughters_Cdk2act,1);
G0durationhist=0;
if G0durationhist
    figure,hist(reltime_Cdk2act/framesperhr,100);
    title('G0 durations for cells entering G1');
end
timesincerisegraph=0;
if timesincerisegraph
    figure, hold on;
    binstep=1; binmin=minlength/framesperhr; binmax=20;
    binscattershade(timesince_Cdk2act/framesperhr,IF_Cdk2act,binstep,binmin,binmax,'b');
    %scatterandfit(timesince_Cdk2act/framesperhr,IF_Cdk2act,'b','r',1,binmin);
    title([IFstring,' levels vs time since Cdk2 activation']);
    xlabel('Time since Cdk2 activation (hrs)'); ylabel([IFstring,' (log2(medianRFU))']);
    xlim([binmin binmax]);
    set(gcf,'color','w','PaperPosition',[0 0 3 3]); %saveas(gcf,'h:\Downloads\Fig3.jpg');
end
actvsnonactgraph=0;
if actvsnonactgraph
    figure, hold on;
    data=[IF_Cdk2nonact;IF_Cdk2act];
    offset=[ones(numel(IF_Cdk2nonact),1);2*ones(numel(IF_Cdk2act),1)];
    boxplot(data,offset,'labels',{'Cdk2 activated','Cdk2 not activated'});
    title([IFstring,' levels with or without Cdk2 activation']);
    ylabel([IFstring,' (log2(medianRFU))']);
    set(gcf,'color','w','PaperPosition',[0 0 3 3]); %saveas(gcf,'h:\Downloads\Fig4.jpg');
end

daughters=~isnan(allstats(:,4));
%%% measure minimum CDK2 activity %%%%%%
daughtertracesCdk2=alltraces1(daughters,:);
daughterstats=allstats(daughters,:);
numtraces=size(daughtertracesCdk2,1);
minval=ones(numtraces,1)*NaN;
maxval=ones(numtraces,1)*NaN;
for i=1:numtraces
    %minval(i)=min(smooth(alltraces1(i,allstats(i,1)+4:allstats(i,1)+minlength-1))); %default 19
    %minval(i)=min(daughtertracesCdk2(i,daughterstats(i,1)+minlength-5:daughterstats(i,1)+minlength-1)); %default 19
    minval(i)=daughtertracesCdk2(i,daughterstats(i,1)+minlength-1);
    maxval(i)=max(daughtertracesCdk2(i,daughterstats(i,1)+minlength-1:daughterstats(i,2)));
end
%figure,hist(maxval,50);
mincutoff=0.5;
maxcutoff=0.6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reltime_mitosis=allstats(daughters,2)-allstats(daughters,1);
reltime_IF=IFdata(daughters);
Cdk2inc=minval>mincutoff & maxval>maxcutoff+0.4;
Cdk2low=minval<mincutoff & maxval<maxcutoff;
reltime_Cdk2inc=reltime_mitosis(Cdk2inc);
reltime_Cdk2low=reltime_mitosis(Cdk2low);
IF_Cdk2inc=reltime_IF(Cdk2inc);
IF_Cdk2low=reltime_IF(Cdk2low);
traces_Cdk2inc=daughtertracesCdk2(Cdk2inc,:);
traces_Cdk2low=daughtertracesCdk2(Cdk2low,:);
incvslowgraph=1;
if incvslowgraph
    figure, hold on;
    for i=1:size(traces_Cdk2inc,1)
        truestart=find(~isnan(traces_Cdk2inc(i,:)),1);
        trueend=find(~isnan(traces_Cdk2inc(i,:)),1,'last');
        line((0:reltime_Cdk2inc(i))/framesperhr,traces_Cdk2inc(i,truestart:trueend),'color','b');
    end
    for i=1:size(traces_Cdk2low,1)
        truestart=find(~isnan(traces_Cdk2low(i,:)),1);
        trueend=find(~isnan(traces_Cdk2low(i,:)),1,'last');
        line((0:reltime_Cdk2low(i))/framesperhr,traces_Cdk2low(i,truestart:trueend),'color','r');
    end
    figure, hold on;
    binstep=1; binmin=minlength/framesperhr; binmax=20;
    scatter(reltime_Cdk2inc/framesperhr,IF_Cdk2inc,'b','o');
    scatter(reltime_Cdk2low/framesperhr,IF_Cdk2low,'r','o','markerfacecolor','r');
    legend(char({'Cdk2-inc';'Cdk2-low'}),'location','southeast');
    %binscattershade(reltime_Cdk2inc/framesperhr,IF_Cdk2inc,binstep,binmin,binmax,'b');
    %binscattershade(reltime_Cdk2low/framesperhr,IF_Cdk2low,binstep,binmin,binmax,'r');
    xlim([minlength/framesperhr binmax]);
    xlabel('Time since Mitosis (hrs)'); ylabel([IFstring,' (log2(medianRFU))']);
    title([IFstring,' levels for Cdk2-inc vs Cdk2-low cells']);
    set(gcf,'color','w','PaperPosition',[0 0 4 4]); saveas(gcf,'h:\Downloads\Fig5.jpg');
    figure, hold on;
    earlyinc=reltime_Cdk2inc/framesperhr<12; IF_Cdk2inc=IF_Cdk2inc(earlyinc);
    earlylow=reltime_Cdk2low/framesperhr<12; IF_Cdk2low=IF_Cdk2low(earlylow);
    data=[IF_Cdk2inc;IF_Cdk2low];
    offset=[ones(numel(IF_Cdk2inc),1);2*ones(numel(IF_Cdk2low),1)];
    boxplot(data,offset,'labels',{'Cdk2-inc','Cdk2-low'});
    title([IFstring,' levels for Cdk2-inc vs Cdk2-low cells']);
    ylabel([IFstring,' (log2(medianRFU))']);
    set(gcf,'color','w','PaperPosition',[0 0 4 4]); saveas(gcf,'h:\Downloads\Fig6.jpg');
end
end


function [tracesCdk2,traces2,IF1,IF2,stats]=main(datadir,shot,minlength)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
load([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
%%% gate traces based on signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(IFdata(:,1))); %take all cells (not just daughters)
%samplecells=find(~isnan(IFdata(:,1)) & ~isnan(genealogy)); %take only daughters
signals=tracedata(samplecells,:,5);
signals2=tracedata(samplecells,:,7)./tracedata(samplecells,:,5);
stats=tracestats(samplecells,:);
minthresh=75; %Sabrina and me:75. Kyuho:150
[samplecells,samplestats]=remove_short_neg_noisy_Cdk2(samplecells,signals,signals2,stats,minlength,minthresh);
%%% gate on another signal (optional) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signals=tracedata(samplecells,:,5);
% %hist(signals(:),100);
% minthresh=0;
% [samplecells,samplestats]=remove_short_neg_noisy(samplecells,signals,samplestats,minlength,minthresh);
%%% calculate signal of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesCdk2=tracedata(:,:,7)./tracedata(:,:,5); %Sabrina:6/5. Kyuho:7/5. HeeWon:7/5.
traces2=tracesCdk2;
stats=tracestats(samplecells,:);
%%% smooth data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(samplecells)
    tracenan=isnan(tracesCdk2(samplecells(i),:));
    tracesCdk2(samplecells(i),:)=smooth(tracesCdk2(samplecells(i),:));
    traces2(samplecells(i),:)=smooth(traces2(samplecells(i),:));
    tracesCdk2(samplecells(i),tracenan)=NaN;
    traces2(samplecells(i),tracenan)=NaN;
end
%%% remove non-sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesCdk2=tracesCdk2(samplecells,:);
traces2=traces2(samplecells,:);
%%% calculate IF values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF1=IFdata(samplecells,8); %Kyuho: pS6nuc=8, pS6ring=9.  Sabrina: p21=7
IF2=IFdata(samplecells,7);
end
