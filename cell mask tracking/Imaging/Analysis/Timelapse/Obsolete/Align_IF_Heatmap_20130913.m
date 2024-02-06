function Align_IF_Heatmap(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces=[];
    allEdU=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [traces,EdU]=main(datadir,shot);
                alltraces=[alltraces;traces];
                allEdU=[allEdU;EdU];
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
[degstarts,degends,badtraces]=getdegstartandend_IF(sampletraces);
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
firstframe=prctile(minframe,10);
latetraces=minframe>=firstframe;
minframe=minframe(latetraces);
sampletraces=sampletraces(latetraces,:);
degstarts=degstarts(latetraces);
degends=degends(latetraces);
minmap=minmap(latetraces,:);
sampleid=sampleid(latetraces);
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[minframe,idx]=sort(minframe);
degstarts=degstarts(idx);
degends=degends(idx);
sampletraces=sampletraces(idx,:);
minmap=minmap(idx,:);
sampleid=sampleid(idx);
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces,1);
sampletraces(minmap==1)=0.5;
for i=1:numtraces
    sampletraces(i,minframe(i):degstarts(i)-1)=1;
    if ~isnan(degends(i))
        sampletraces(i,degstarts(i):degends(i)-1)=0;
        sampletraces(i,degends(i):end)=1;
    else
        sampletraces(i,degstarts(i):end)=0;
    end
end
traceid=1:numtraces;
time=((firstframe:numframes)-numframes)/framesperhr;
data=sampletraces(:,firstframe:numframes);
imagesc(time,traceid,data);
bluecode=[0 0 1]; blackcode=[0 0 0]; yellowcode=[1 1 0];
cmap=makecmap(bluecode,blackcode,yellowcode);
colormap(cmap);
%%% overlay markers at timepoints of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hold on;
% scatter((degstarts-numframes)/framesperhr,traceid,4,'r','*');
% notnan=~isnan(degends);
% scatter((degends(notnan)-numframes)/framesperhr,traceid(notnan),4,'g','*');
%%% overlay EdU values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EdUvals=log(allEdU(sampleid));
maxEdU=prctile(EdUvals,90);
minEdU=prctile(EdUvals,10);
EdUnorm=(EdUvals-minEdU)/(maxEdU-minEdU);
EdUnorm(EdUnorm>1)=1;
EdUnorm(EdUnorm<0)=0;
EdUnorm=round(EdUnorm*128);
EdUnorm(EdUnorm==0)=1;
EdUtime=ones(numtraces,1)*-0.2;
redcode=[1 0 0]; blackcode=[0 0 0]; greencode=[0 1 0];
EdUcmap=makecmap(redcode,blackcode,greencode);
EdUcolor=EdUcmap(EdUnorm,:);
scatter(EdUtime,traceid,20,EdUcolor,'*');
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Time relative to Immunostain (hrs)'); ylabel('Individual traces (ordered by mitosis)');
set(gcf,'color','w','PaperPosition',[0 0 5 8]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%%% IF by category %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ligaseoff=find(sampletraces(:,end)==1);
ligaseon=find(sampletraces(:,end)==0);
edudata=[EdUvals(ligaseoff);EdUvals(ligaseon)];
offset=[ones(length(ligaseoff),1);2*ones(length(ligaseon),1)];
figure,boxplot(axes,edudata,offset,'labels',{'sensor-high','sensor-low'});
set(gca,'Fontsize',14);
h=findobj(gca,'Type','text');
set(h,'FontSize',14)
set(gcf,'color','w');
end


function [traces,EdUval]=main(datadir,shot)
sensorchannel=6;
EdUchannel=10; %10=max
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
load([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
totalframes=size(tracedata,2);
%%% find daughters that are stained %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(IFdata(:,1)) & ~isnan(tracestats(:,2)));
samplesignals=tracedata(samplecells,:,sensorchannel);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
[samplecells,samplesignals]=gatetraces(samplecells,samplesignals,tracestats);
%%% align traces to stain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces=ones(numel(samplecells),totalframes)*NaN;
for i=1:numel(samplecells)
    s=samplecells(i);
    frames=tracestats(s,1):tracestats(s,3);
    numframes=tracestats(s,3)-tracestats(s,1)+1;
    traces(i,end-numframes+1:end)=samplesignals(i,frames);
end
EdUval=IFdata(samplecells,EdUchannel);
end
