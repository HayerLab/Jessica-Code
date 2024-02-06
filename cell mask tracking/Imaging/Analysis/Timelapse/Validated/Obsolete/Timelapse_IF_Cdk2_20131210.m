function Timelapse_IF_Cdk2(conditions,datadir)
condnum=size(conditions,1);
minlength=40;
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
maxlength=prctile(allstats(:,3),90);
sampleid=(1:numtraces)';
detectrisetime=1;
if detectrisetime
    %%% detect onset of CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relative=0; %returns Cdk2start time relative to trace start
    [~,Cdk2start,~,badtraces]=getCdk2features_mother(alltraces1,allstats,minlength,relative);
    Cdk2start(isnan(Cdk2start))=0;
    sampleid(badtraces>0)=[];
    %%% gate by trace length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tracelength=allstats(:,3);
    maxlength=prctile(tracelength,100);
    longtraces=tracelength>maxlength;
    sampleid(longtraces)=[];
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
    minval(i)=min(smooth(alltraces1(i,allstats(i,1)+4:allstats(i,1)+19)));
end
%%% remove IF outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=allIF1;
IFdata(IFdata<1)=1;
IFdata=log2(IFdata);
%hist(IFdata(:),100);
[mu1,sig1]=normfit(IFdata);
outliers=IFdata>(mu1+3*sig1) | IFdata<(mu1-3*sig1); % | IFdata<1
IFdata(outliers)=[];
alltraces1(outliers,:)=[];
alltraces2(outliers,:)=[];
Cdk2start(outliers)=[];
allstats(outliers,:)=[];
%%% align data to immunostain time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltraces1,1);
alignedtraces1=ones(numtraces,numframes)*NaN;
alignedtraces2=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    alignedtraces1(i,numframes-allstats(i,3)+1:numframes)=alltraces1(i,allstats(i,1):allstats(i,2));
    alignedtraces2(i,numframes-allstats(i,3)+1:numframes)=alltraces2(i,allstats(i,1):allstats(i,2));
end
alignedtraces1=alignedtraces1(:,numframes-maxlength+1:numframes);
alignedtraces2=alignedtraces2(:,numframes-maxlength+1:numframes);
%%% display graphs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracedata=alignedtraces1;
time=(numframes-maxlength-1:numframes)/framesperhr;
heatmapoption=1;
if heatmapoption
    tracedata(tracedata>2.5)=2.5;
    tracedata(tracedata<0.3)=0.3;
    tracedata(isnan(tracedata))=0.3; %make invalid data black
    traceid=1:numtraces;
    imagesc(time,traceid,tracedata);
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
        onlyreal=find(Cdk2start~=0);
        Cdk2start_onlyreal=Cdk2start(onlyreal);
        traceid_onlyreal=traceid(onlyreal);
        scatter((Cdk2start_onlyreal)/framesperhr,traceid_onlyreal,8,'r','*');
    end
    %%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xlabel('Time relative to Mitosis (hrs)'); ylabel('Individual traces');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]);
    saveas(gcf,'h:\Downloads\Fig.jpg');
    %%% make IF panel to accompany heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure,imagesc(1,traceid,IFdata);
    set(gca,'XTick',[],'YTick',[]);
    set(gcf,'color','w','PaperPosition',[0 0 1 8]);
    saveas(gcf,'h:\Downloads\Fig2.jpg');
end

%%% G1 timelapse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keyboard;
figure, hold on;
scatterandfit(G1time(G1inc),G1IF(G1inc),'b','b',1,4);
scatterandfit(G1time(G1low),G1IF(G1low),'r','r',1,4);
title('Cyclin D1 levels in G1: Cdk2-inc vs Cdk2-low'); xlabel('Time since mitosis (hrs)'); ylabel('Cyclin D1 (log(median-RFU))');
axis([3 15 3 8]);
set(gcf,'color','w','PaperPosition',[0 0 6 3]); saveas(gcf,'h:\Downloads\Fig1.jpg');
%%% G1 boxplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1data=[G1IF(G1inc);G1IF(G1low)];
G1offset=[ones(sum(G1inc),1);2*ones(sum(G1low),1)];
figure,boxplot(axes,G1data,G1offset,'labels',{'Cdk2-inc','Cdk2-low'});
title('Cyclin D1 levels in G1 phase'); xlabel('Cdk2-inc vs Cdk2-low)'); ylabel('Cyclin D1 (log(median-RFU))');
set(gcf,'color','w','PaperPosition',[0 0 3 3]); saveas(gcf,'h:\Downloads\Fig2.jpg');
%%% S and G2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scatter(Stime(Sinc),SIF(Sinc),'b.');
scatter(Stime(Slow),SIF(Slow),'r.');
scatter(SG2time(G2inc),G2IF(G2inc),'b+');
scatter(SG2time(G2low),G2IF(G2low),'r+');
axis([0 20 2 8]);
title('Cyclin D1 levels: Cdk2-inc vs Cdk2-low'); xlabel('Time since G1/S (hrs)'); ylabel('Cyclin D1 (log(median-RFU))');
set(gcf,'color','w','PaperPosition',[0 0 8 3]); saveas(gcf,'h:\Downloads\Fig1.jpg');
%%% box plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sdata=[SIF(Sinc);SIF(Slow)];
Soffset=[ones(sum(Sinc),1);2*ones(sum(Slow),1)];
figure,boxplot(axes,Sdata,Soffset,'labels',{'Cdk2-inc','Cdk2-low'});
title('Cyclin D1 levels in S phase'); xlabel('Cdk2-inc vs Cdk2-low)'); ylabel('Cyclin D1 (log(median-RFU))');
set(gcf,'color','w','PaperPosition',[0 0 3 3]); saveas(gcf,'h:\Downloads\Fig2.jpg');
G2data=[G2IF(G2inc);G2IF(G2low)];
G2offset=[ones(sum(G2inc),1);2*ones(sum(G2low),1)];
figure,boxplot(axes,G2data,G2offset,'labels',{'Cdk2-inc','Cdk2-low'});
title('Cyclin D1 levels in G2'); xlabel('Cdk2-inc vs Cdk2-low)'); ylabel('Cyclin D1 (log(median-RFU))');
set(gcf,'color','w','PaperPosition',[0 0 3 3]); saveas(gcf,'h:\Downloads\Fig3.jpg');
end


function [traces1,traces2,IF1,IF2,stats]=main(datadir,shot,minlength)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
load([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
%%% gate traces based on signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(IFdata(:,1)));
signals=tracedata(samplecells,:,5);
stats=tracestats(samplecells,:);
%minlength=40;
%hist(signals(:),100);
minthresh=100; %200 for Sabrina and Kyuho
[samplecells,~]=remove_short_neg_noisy(samplecells,signals,stats,minlength,minthresh);
%%% gate on another signal (optional) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signals=tracedata(samplecells,:,5);
% stats=tracestats(samplecells,:);
% minthresh=200;
% [samplecells,samplestats]=remove_short_neg_noisy(samplecells,signals,stats,minlength,minthresh);
%%% calculate signal of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces1=tracedata(samplecells,:,7)./tracedata(samplecells,:,5);
traces2=tracedata(samplecells,:,6);
stats=tracestats(samplecells,:);
%%% calculate IF values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF1=IFdata(samplecells,8); %ring of w4
IF2=IFdata(samplecells,8); %nuc of w5
end
