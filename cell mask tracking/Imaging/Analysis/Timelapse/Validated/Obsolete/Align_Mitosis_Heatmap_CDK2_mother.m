function Align_Mitosis_Heatmap_CDK2(conditions,datadir,motheroption)
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
[~,Cdk2start,~,badtraces]=getCdk2features_mother(alltraces,allstats,relative);
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
data=data(:,1:maxframe);
%%% display mother data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mothermax=0; datamother=[];
if motheroption
    datamother=ones(numtraces,numframes)*NaN;
    for i=1:numtraces
        firstframe=find(~isnan(alltraces(i,:)),1,'first');
        lastframe=allstats(i,1)-1;
        duration=lastframe-firstframe+1;
        datamother(i,numframes-duration+1:numframes)=alltraces(i,firstframe:lastframe);
    end
    mothermax=40;
    datamother=datamother(:,numframes-mothermax+1:numframes);
end
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datatotal=[datamother data];
datatotal(datatotal>2)=2;
datatotal(datatotal<0.2)=0.2;
datatotal(isnan(datatotal))=0.2; %make invalid data black
traceid=1:numtraces;
time=(-mothermax+1:maxframe)/framesperhr;
imagesc(time,traceid,datatotal);
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
    scatter((Cdk2start+mothermax)/framesperhr,traceid,8,'r','*');
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
samplecells=find(~isnan(tracestats(:,4)));
%%% gate traces based on signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals=tracedata(samplecells,:,5);
stats=tracestats(samplecells,:);
minlength=40;
minthresh=200;
[samplecells,samplestats]=remove_short_neg_noisy(samplecells,signals,stats,minlength,minthresh);
%%% calculate signal of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces=tracedata(:,:,7)./tracedata(:,:,5);
%%% append mother data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(samplecells)
    prevframe=samplestats(i,1)-1;
    traces(samplecells(i),1:prevframe)=traces(samplestats(i,4),1:prevframe);
end
%%% remove non-sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces=traces(samplecells,:);
end