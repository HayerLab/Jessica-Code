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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
[~,numframes]=size(alltraces);
sampletraces=alltraces;
%%% detect degradation start and end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[degstarts,degends,badtraces]=getdegstartandend(sampletraces);
sampletraces(badtraces>0,:)=[];
degstarts(badtraces>0)=[];
degends(badtraces>0)=[];
%%% gate by trace length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces,1);
minmap=zeros(numtraces,numframes);
maxframe=zeros(numtraces,1);
for i=1:numtraces
    maxframe(i)=find(~isnan(sampletraces(i,:)),1,'last');
    minmap(i,maxframe(i)+1:end)=1;
end
lastframe=prctile(maxframe,90);
shorttraces=maxframe<=lastframe;
maxframe=maxframe(shorttraces);
sampletraces=sampletraces(shorttraces,:);
degstarts=degstarts(shorttraces);
degends=degends(shorttraces);
minmap=minmap(shorttraces,:);
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sorttraces=1;
if sorttraces
    [degstarts,idx]=sort(degstarts);
    maxframe=maxframe(idx);
    degends=degends(idx);
    sampletraces=sampletraces(idx,:);
    minmap=minmap(idx,:);
end
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces,1);
%%%%%% convert continuous to binary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
continuous2binary=0;
if continuous2binary
    for i=1:numtraces
        endframe=find(~isnan(sampletraces(i,:)),1,'last');
        sampletraces(i,1:degstarts(i)-1)=1;
        if ~isnan(degends(i))
            sampletraces(i,degstarts(i):degends(i)-1)=0;
            sampletraces(i,degends(i):endframe)=1;
        else
            sampletraces(i,degstarts(i):endframe)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampletraces(minmap==1)=0.5;
traceid=1:numtraces;
time=(1:lastframe)/framesperhr;
data=sampletraces(:,1:lastframe);
imagesc(time,traceid,data);
bluecode=[0 0 1]; blackcode=[0 0 0]; yellowcode=[1 1 0];
cmap=makecmap(bluecode,blackcode,yellowcode);
colormap(cmap);
%%% overlay markers at timepoints of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addmarks=1;
if addmarks
    hold on;
    scatter(degstarts/framesperhr,traceid,8,'r','*');
    scatter(degends/framesperhr,traceid,8,'g','*');
end
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Time relative to Mitosis (hrs)'); ylabel('Individual traces (sorted by degradation time)');
set(gcf,'color','w','PaperPosition',[0 0 4 7]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end


function traces=main(datadir,shot)
sensorchannel=6;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalframes=size(tracedata,2);
samplecells=find(~isnan(tracestats(:,2)));
%%% get cells that are also mothers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ismother=zeros(numel(samplecells),1);
for i=1:numel(samplecells)
    s=samplecells(i);
    ismother(i)=ismember(s,tracestats(:,2));
end
samplecells=samplecells(ismother>0);
%%% get data from sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplesignals=tracedata(samplecells,:,sensorchannel);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
[~,samplesignals,samplestats]=gatetraces_pipdegron(samplecells,samplesignals,samplestats);
%%% align traces to mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(samplesignals,1);
traces=ones(numtraces,totalframes)*NaN;
for i=1:numtraces
    frames=samplestats(i,1):samplestats(i,3);
    traces(i,1:numel(frames))=samplesignals(i,frames);
end
end

%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
sample=523;
traces=sampletraces(sample,:);
starts=degstarts(sample);
ends=degends(sample);
maxf=maxframe(sample);
samplesize=size(traces,1);
for i=1:samplesize
    %figure(ceil(i/24)); set(gcf,'color','w');
    %subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    x=1:maxf(i);
    plot(x/framesperhr,traces(i,x),'linewidth',4);
    xlim([1 21.6]);
    hold on;
    plot(starts(i)/framesperhr,traces(i,starts(i)),'ro','markerfacecolor','r','markersize',8);
    plot(ends(i)/framesperhr,traces(i,ends(i)),'go','markerfacecolor','g','markersize',8);
end
xlabel('Time relative to Mitosis (hrs)'); ylabel('sensor level');
set(gcf,'color','w','PaperPosition',[0 0 3 2]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%