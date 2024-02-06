function Timelapse_IF_Cdk2_pipdeg(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    all_dhb=[];
    all_pip=[];
    all_IF=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [traces_dhb,traces_pip,IF]=main(datadir,shot);
                all_dhb=[all_dhb;traces_dhb];
                all_pip=[all_pip;traces_pip];
                all_IF=[all_IF;IF];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
[numtraces,numframes]=size(all_dhb);
sampletraces_dhb=all_dhb;
sampletraces_pip=all_pip;
sampleid=(1:numtraces)';
%%% extract Cdk2 features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstframe=ones(numtraces,1)*NaN;
minval=ones(numtraces,1)*NaN;
for i=1:numtraces
    firstframe(i)=find(~isnan(sampletraces_dhb(i,:)),1);
    minval(i)=sampletraces_dhb(i,firstframe(i)+19);
end
badtraces=minval<0.2 | minval>1.5;
firstframe(badtraces)=[];
minval(badtraces)=[];
sampletraces_pip(badtraces,:)=[];
sampleid(badtraces)=[];
%%% extract pipdeg features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relative=0;
[degstarts,degends,badtraces]=getdegstartandend_openended(sampletraces_pip,relative);
degstarts(badtraces)=[];
degends(badtraces)=[];
firstframe(badtraces)=[];
minval(badtraces)=[];
sampleid(badtraces)=[];
%%% remove IF outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleIF=all_IF(sampleid);
[mu,sig]=normfit(sampleIF);
outliers=sampleIF>(mu+3*sig) | sampleIF<(mu-3*sig) | sampleIF<1;
sampleIF=sampleIF(~outliers);
degstarts=degstarts(~outliers);
degends=degends(~outliers);
firstframe=firstframe(~outliers);
minval=minval(~outliers);
EdUvals=log(sampleIF);
%%% Calc G1-phase IF vs time by dhb ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1sample=isnan(degstarts);
G1start=firstframe(G1sample);
G1IF=EdUvals(G1sample);
G1time=(numframes-G1start)/framesperhr;
G1minval=minval(G1sample);
G1inc=G1minval>0.56;
G1low=G1minval<=0.56;
%%% Calc S-phase IF vs time by dhb ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ssample=~isnan(degstarts) & isnan(degends);
Sdegstarts=degstarts(Ssample);
SIF=EdUvals(Ssample);
Stime=(numframes-Sdegstarts)/framesperhr;
Sminval=minval(Ssample);
Sinc=Sminval>0.56;
Slow=Sminval<=0.56;
%%% Calc G2-phase IF vs time by dhb ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G2sample=~isnan(degstarts) & ~isnan(degends);
G2degstarts=degstarts(G2sample);
G2degends=degends(G2sample);
G2time=(numframes-G2degends)/framesperhr;
SG2time=(numframes-G2degstarts)/framesperhr;
G2IF=EdUvals(G2sample);
G2minval=minval(G2sample);
G2inc=G2minval>0.56;
G2low=G2minval<=0.56;
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
%%% G2 timelapse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, hold on;
scatterandfit(G2time(G2inc),G2IF(G2inc),'b','b',1,1);
scatterandfit(G2time(G2low),G2IF(G2low),'r','r',1,1);
title('Cyclin D1 levels in G2: Cdk2-inc vs Cdk2-low'); xlabel('Time since S/G2 (hrs)'); ylabel('Cyclin D1 (log(median-RFU))');
axis([0 2.5 4 7.5]);
set(gcf,'color','w','PaperPosition',[0 0 3 3]); saveas(gcf,'h:\Downloads\Fig1.jpg');
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


function [traces_dhb,traces_pip,IFval]=main(datadir,shot)
pipchannel=6;
IFchannel=9; %9=med; 10=max 11=mean
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
load([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
totalframes=size(tracedata,2);
%%% find daughters that are stained %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(IFdata(:,1)) & ~isnan(tracestats(:,2)));
samplesignals=tracedata(samplecells,:,5);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
samplecells=gatetraces_cdk2(samplecells,samplesignals,samplestats,20);
samplesignals=tracedata(samplecells,:,pipchannel);
samplestats=tracestats(samplecells,:);
[samplecells,sig_pip,~]=gatetraces_pipdegron(samplecells,samplesignals,samplestats);
%%% get data from sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig_dhb=tracedata(samplecells,:,7)./tracedata(samplecells,:,5);
sig_dhb(sig_dhb<0)=0;
%%% align traces to stain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces_pip=ones(numel(samplecells),totalframes)*NaN;
traces_dhb=ones(numel(samplecells),totalframes)*NaN;
for i=1:numel(samplecells)
    s=samplecells(i);
    frames=tracestats(s,1):tracestats(s,3);
    traces_pip(i,end-numel(frames)+1:end)=sig_pip(i,frames);
    traces_dhb(i,end-numel(frames)+1:end)=sig_dhb(i,frames);
end
IFval=IFdata(samplecells,IFchannel);
end
