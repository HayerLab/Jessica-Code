function Align_IF_Heatmap(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces=[];
    allIF=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [traces,IF]=main(datadir,shot);
                alltraces=[alltraces;traces];
                allIF=[allIF;IF];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
[numtraces,numframes]=size(alltraces);
sampletraces=alltraces;
sampleid=(1:numtraces)';
%%% detect degradation start and end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[degstarts,degends,badtraces]=getdegstartandend(sampletraces);
sampletraces(badtraces>0,:)=[];
degstarts(badtraces>0)=[];
degends(badtraces>0)=[];
sampleid(badtraces>0)=[];
%%% gate by trace length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces,1);
minmap=zeros(numtraces,numframes);
minframe=zeros(numtraces,1);
for i=1:numtraces
    minframe(i)=find(~isnan(sampletraces(i,:)),1,'first');
    minmap(i,1:minframe(i)-1)=1;
end
firstframe=prctile(minframe,0);
latetraces=minframe>=firstframe;
minframe=minframe(latetraces);
sampletraces=sampletraces(latetraces,:);
degstarts=degstarts(latetraces);
degends=degends(latetraces);
minmap=minmap(latetraces,:);
sampleid=sampleid(latetraces);
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sorttraces=0;
if sorttraces
    [minframe,idx]=sort(minframe);
    degstarts=degstarts(idx);
    degends=degends(idx);
    sampletraces=sampletraces(idx,:);
    minmap=minmap(idx,:);
    sampleid=sampleid(idx);
end
%%% remove IF outliers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleIF=allIF(sampleid);
[mu,sig]=normfit(sampleIF);
outliers=sampleIF>mu+3*sig | sampleIF<mu-3*sig | sampleIF<1;
sampleIF=sampleIF(~outliers);
degstarts=degstarts(~outliers);
degends=degends(~outliers);
sampletraces=sampletraces(~outliers,:);
minmap=minmap(~outliers,:);
sampleid=sampleid(~outliers);
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces,1);
%%%%%% convert continuous to binary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
continuous2binary=0;
if continuous2binary
    for i=1:numtraces
        sampletraces(i,minframe(i):degstarts(i)-1)=1;
        if ~isnan(degends(i))
            sampletraces(i,degstarts(i):degends(i)-1)=0;
            sampletraces(i,degends(i):end)=1;
        else
            sampletraces(i,degstarts(i):end)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampletraces(minmap==1)=0.5;
traceid=1:numtraces;
time=((firstframe:numframes)-numframes)/framesperhr;
data=sampletraces(:,firstframe:numframes);
imagesc(time,traceid,data);
bluecode=[0 0 1]; blackcode=[0 0 0]; yellowcode=[1 1 0];
cmap=makecmap(bluecode,blackcode,yellowcode);
colormap(cmap);
%%% overlay markers at timepoints of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addmarks=0;
if addmarks
    hold on;
    scatter((degstarts-numframes)/framesperhr,traceid,4,'r','*');
    notnan=~isnan(degends);
    scatter((degends(notnan)-numframes)/framesperhr,traceid(notnan),4,'g','*');
end
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Time relative to Immunostain (hrs)'); ylabel('Individual traces (ordered by mitosis)');
set(gcf,'color','w','PaperPosition',[0 0 5 8]);
saveas(gcf,'h:\Downloads\Fig1.jpg');
%%% overlay EdU values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sampleIF(sampleIF<1)=1;
EdUvals=log(sampleIF);
% maxEdU=prctile(EdUvals,90);
% minEdU=prctile(EdUvals,10);
% EdUnorm=(EdUvals-minEdU)/(maxEdU-minEdU);
% EdUnorm(EdUnorm>1)=1;
% EdUnorm(EdUnorm<0)=0;
EdUnorm=EdUvals;
%%% make EdU panel to accompany heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,imagesc(1,traceid,EdUnorm);
set(gca,'XTick',[],'YTick',[]);
set(gcf,'color','w','PaperPosition',[0 0 1 8]);
saveas(gcf,'h:\Downloads\Fig2.jpg');
%%% IF by category %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ligaseoff=find(sampletraces(:,end)==1);
ligaseon=find(sampletraces(:,end)==0);
edudata=[EdUvals(ligaseoff);EdUvals(ligaseon)];
offset=[ones(length(ligaseoff),1);2*ones(length(ligaseon),1)];
figure,boxplot(axes,edudata,offset,'labels',{'sensor-high','sensor-low'});
set(gca,'Fontsize',14);
h=findobj(gca,'Type','text');
set(h,'FontSize',14)
set(gcf,'color','w','PaperPosition',[0 0 5 8]);
saveas(gcf,'h:\Downloads\Fig3.jpg');
%%% plot by time since S and G2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ssample=~isnan(degstarts) & isnan(degends);
Sdegstarts=degstarts(Ssample);
SIF=EdUnorm(Ssample);
Stime=(numframes-Sdegstarts)/framesperhr;
scatter(Stime,SIF,'ro');
G2sample=~isnan(degstarts) & ~isnan(degends);
G2degstarts=degstarts(G2sample);
G2IF=EdUnorm(G2sample);
SG2time=(numframes-G2degstarts)/framesperhr;
hold on; scatter(SG2time,G2IF,'bo');
keyboard;
%%% overlay mean and errorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SG2time=[Stime;SG2time];
SG2IF=[SIF;G2IF];
maxbin=ceil(max(SG2time));
if mod(maxbin,2)==0
    maxbin=maxbin-1;
end
for i=1:2:maxbin
    binsample=SG2IF(SG2time>i-1 & SG2time<=i+1);
    if isempty(binsample)
        continue;
    end
    [errmean,errstd]=normfit(binsample);
    errorbar(i,errmean,errstd,'-k','linewidth',2);
end
keyboard;
axis([0 25 3 10]);
xlabel('Time after G1/S (hrs)'); ylabel('EdU Incorporation (log(median-RFU))');
set(gcf,'color','w','PaperPosition',[0 0 10 5]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%%% plot by time since S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ssample=~isnan(degstarts) & isnan(degends);
Sdegstarts=degstarts(Ssample);
SIF=EdUnorm(Ssample);
Stime=(numframes-Sdegstarts)/framesperhr;
scatter(Stime,SIF,'ro');
axis([0 12 2.5 7]);
%%% plot by time since G2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G2sample=~isnan(degends);
G2degends=degends(G2sample);
G2IF=EdUnorm(G2sample);
G2time=numframes-G2degends;
figure,scatter(G2time/framesperhr,G2IF);
axis([0 10 4 8]);
end


function [traces,IFval]=main(datadir,shot)
sensorchannel=6;
IFchannel=9; %9=med; 10=max 11=mean
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
load([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
totalframes=size(tracedata,2);
%%% find daughters that are stained %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(IFdata(:,1)) & ~isnan(tracestats(:,2)));
samplesignals=tracedata(samplecells,:,sensorchannel);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
[samplecells,samplesignals]=gatetraces_pipdegron(samplecells,samplesignals,samplestats);
%%% align traces to stain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces=ones(numel(samplecells),totalframes)*NaN;
for i=1:numel(samplecells)
    s=samplecells(i);
    frames=tracestats(s,1):tracestats(s,3);
    traces(i,end-numel(frames)+1:end)=samplesignals(i,frames);
end
IFval=IFdata(samplecells,IFchannel);
end
