function Align_Mitosis_Heatmap_CDK2(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces=[];
    allstats=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                %shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [traces,stats]=main(datadir,shot);
                alltraces=[alltraces;traces];
                allstats=[allstats;stats];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
numframes=size(alltraces,2);
%%% detect onset of CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relative=1; %returns Cdk2start time relative to mitosis
[~,Cdk2start,~,badtraces]=getCdk2features(alltraces,relative);
Cdk2start(isnan(Cdk2start))=0;
Cdk2start(badtraces>0)=[];
alltraces(badtraces>0,:)=[];
allstats(badtraces>0,:)=[];
%%% gate by trace length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracelength=allstats(:,3);
maxframe=prctile(tracelength,90);
shorttraces=tracelength<=maxframe;
Cdk2start=Cdk2start(shorttraces);
alltraces=alltraces(shorttraces,:);
allstats=allstats(shorttraces,:);
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sorttraces=1;
if sorttraces
    [Cdk2start,idx]=sort(Cdk2start);
    alltraces=alltraces(idx,:);
    allstats=allstats(idx,:);
end
%%% align data to mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltraces,1);
data=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    data(i,1:allstats(i,3))=alltraces(i,allstats(i,1):allstats(i,2));
end
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=data(:,1:maxframe);
data(isnan(data))=0.2; %make invalid data black
traceid=1:numtraces;
time=(1:maxframe)/framesperhr;
imagesc(time,traceid,data);
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
    scatter(Cdk2start/framesperhr,traceid,8,'r','*');
end
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Time relative to Mitosis (hrs)'); ylabel('Individual traces');
set(gcf,'color','w','PaperPosition',[0 0 4 7]);
saveas(gcf,'h:\Downloads\Fig.jpg');
fprintf('%0.0f traces\n',numtraces);
end


function [traces,samplestats]=main(datadir,shot)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daughtercells=find(~isnan(tracestats(:,2)));
%%% gate traces based on signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals=tracedata(daughtercells,:,5);
stats=tracestats(daughtercells,:);
minlength=40;
minthresh=200;
[samplecells,samplestats]=remove_short_neg_noisy(daughtercells,signals,stats,minlength,minthresh);
%%% get data from sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces=tracedata(samplecells,:,7)./tracedata(samplecells,:,5);
traces(traces<0.2)=0.2;
traces(traces>2)=2;
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