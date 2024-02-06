function Align_Mitosis_Heatmap(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                traces=main(datadir,shot);
                alltraces=[alltraces;traces];
            end
        end
    end
end
%%% normalize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
[numtraces,numframes]=size(alltraces);
maxframe=zeros(numtraces,1);
alltracesnorm=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    maxframe(i)=find(~isnan(alltraces(i,:)),1,'last');
    maxval=max(alltraces(i,:));
    alltracesnorm(i,:)=alltraces(i,:)/maxval; %normalize
end
%%% gate by trace length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lastframe=prctile(maxframe,50);
longtraces=maxframe>=lastframe;
sampletraces=alltracesnorm(longtraces,1:lastframe);
%%% detect degradation start and end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[degstarts,degends,badtraces]=getdegstartandend(sampletraces);
sampletraces(badtraces>0,:)=[];
degstarts(badtraces>0)=[];
degends(badtraces>0)=[];
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,idx]=sort(degstarts);
degstarts=degstarts(idx);
degends=degends(idx);
sampletraces=sampletraces(idx,:);
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces,1);
traceid=1:numtraces;
imagesc((1:lastframe)/framesperhr,traceid,sampletraces);
bluecode=[0 0 1]; blackcode=[0 0 0]; yellowcode=[1 1 0];
cmap=makecmap(bluecode,blackcode,yellowcode);
colormap(cmap);
%%% overlay markers at timepoints of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
scatter(degstarts/framesperhr,traceid,4,'r','*');
notnan=~isnan(degends);
scatter(degends(notnan)/framesperhr,traceid(notnan),4,'g','*');
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Time relative to mitosis (hrs)'); ylabel('Individual traces (ordered by degradation time)');
set(gcf,'color','w','PaperPosition',[0 0 5 8]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end


function traces=main(datadir,shot)
sensorchannel=6;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are mothers and daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalframes=size(tracedata,2);
samplecells=find(~isnan(tracestats(:,2)));
ismother=zeros(numel(samplecells),1);
for i=1:numel(samplecells)
    s=samplecells(i);
    ismother(i)=ismember(s,tracestats(:,2));
end
samplecells=samplecells(ismother>0);
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numsamples=numel(samplecells);
noise=max(diff(tracedata(samplecells,:,sensorchannel),2,2),[],2)>150;
firstandlast=zeros(numsamples,1);
for i=1:numsamples
    s=samplecells(i);
    firstandlast(i)=tracedata(s,tracestats(s,1),sensorchannel)>100 ...
        && tracedata(s,tracestats(s,3),sensorchannel)>100 ...
        && min(tracedata(s,:,sensorchannel))<50;
end
samplecells=samplecells(~noise & firstandlast);
%%% align to mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces=ones(numel(samplecells),totalframes)*NaN;
for i=1:numel(samplecells)
    s=samplecells(i);
    numframes=tracestats(s,3)-tracestats(s,1)+1;
    traces(i,1:numframes)=tracedata(s,tracestats(s,1):tracestats(s,3),sensorchannel);
end
end
