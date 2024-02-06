function Align_Mitosis_Heatmap_CDK2(conditions,datadir,traceoption)
motheroption=1;
minlength=20;
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltracesCdk2=[];
    alltraces2=[];
    allstats=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [tracesCdk2,traces2,stats]=main(datadir,shot,minlength);
                alltracesCdk2=[alltracesCdk2;tracesCdk2];
                alltraces2=[alltraces2;traces2];
                allstats=[allstats;stats];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
[numtraces,numframes]=size(alltracesCdk2);
sampleid=(1:numtraces)';
detectrisetime=1;
if detectrisetime
    %%% detect onset of CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relative=1; %returns Cdk2start time relative to mitosis
    [~,Cdk2start,~,badtraces]=getCdk2features_mother(alltracesCdk2,allstats,minlength,relative);
    Cdk2start(isnan(Cdk2start))=0;
    sampleid(badtraces>0)=[];
end
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sorttraces=1;
if sorttraces
    [Cdk2start,idx]=sort(Cdk2start(sampleid));
    ordered=sampleid(idx);
    alltracesCdk2=alltracesCdk2(ordered,:);
    alltraces2=alltraces2(ordered,:);
    allstats=allstats(ordered,:);
end
%%% measure minimum CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltracesCdk2,1);
minval=ones(numtraces,1)*NaN;
maxval=ones(numtraces,1)*NaN;
for i=1:numtraces
    %minval(i)=min(smooth(alltracesCdk2(i,allstats(i,1)+4:allstats(i,1)+minlength-1))); %default 19
    %minval(i)=min(alltracesCdk2(i,allstats(i,1)+4:allstats(i,1)+minlength-1)); %default 19
    minval(i)=smooth(alltracesCdk2(i,allstats(i,1)+minlength-1)); %Sabrina uses 2hr
    maxval(i)=max(smooth(alltracesCdk2(i,allstats(i,1)+minlength-1:allstats(i,2)))); %max from minlength-end
end

%figure,hist(maxval,50);
mincutoff=0.65; %minval at 20frames
maxcutoff=1.0; %maxval cutoff F20-end
%%% align data to mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltracesCdk2,1);
dataCdk2=ones(numtraces,numframes)*NaN;
data2=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    dataCdk2(i,1:allstats(i,3))=alltracesCdk2(i,allstats(i,1):allstats(i,2));
    data2(i,1:allstats(i,3))=alltraces2(i,allstats(i,1):allstats(i,2));
end
maxlength=prctile(allstats(:,3),90);
dataCdk2=dataCdk2(:,1:maxlength);
data2=data2(:,1:maxlength);
%%% align mother data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mothermax=0; dataCdk2mother=[]; data2mother=[];
if motheroption
    dataCdk2mother=ones(numtraces,numframes)*NaN;
    data2mother=ones(numtraces,numframes)*NaN;
    for i=1:numtraces
        firstframe=find(~isnan(alltracesCdk2(i,:)),1,'first');
        lastframe=allstats(i,1)-1;
        duration=lastframe-firstframe+1;
        dataCdk2mother(i,numframes-duration+1:numframes)=alltracesCdk2(i,firstframe:lastframe);
        data2mother(i,numframes-duration+1:numframes)=alltraces2(i,firstframe:lastframe);
    end
    mothermax=80;
    dataCdk2mother=dataCdk2mother(:,numframes-mothermax+1:numframes);
    data2mother=data2mother(:,numframes-mothermax+1:numframes);
end
%%% concatenate mother & daughter data (mix & match if desired) %%%%%%%%%%%
if traceoption==1
    datatotal=[dataCdk2mother dataCdk2];
elseif traceoption==2
    datatotal=[data2mother data2];
end
%datatotal=[data2mother dataCdk2];
time=(-mothermax+1:maxlength)/framesperhr;
%%% generate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cdk2inc=minval>mincutoff & maxval>maxcutoff;
Cdk2low=minval<mincutoff & maxval<maxcutoff;
dataCdk2inc=datatotal(Cdk2inc,:);
dataCdk2low=datatotal(Cdk2low,:);
durationinc=allstats(Cdk2inc,3)/framesperhr;
durationlow=allstats(Cdk2low,3)/framesperhr;
%figure,hist(durationinc,2:1:24);saveas(gcf,'h:\Downloads\Fig1.jpg');
%figure,hist(durationlow,2:1:24);saveas(gcf,'h:\Downloads\Fig2.jpg');
traceoption=1;
if traceoption
    figure, hold on;
    for i=1:size(dataCdk2inc,1)
        line(time,dataCdk2inc(i,:),'color','b');
    end
    for i=1:size(dataCdk2low,1)
        line(time,dataCdk2low(i,:),'color','r');
    end
    xlabel('Time relative to Mitosis (hrs)');
    ylabel('CDK2 activity'); ylabel('ERK activity');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]);
    saveas(gcf,'h:\Downloads\Fig1.jpg');
    
    figure, hold on;
    binstep=1;
    bincurveshade(time,dataCdk2inc,binstep,'b');
    bincurveshade(time,dataCdk2low,binstep,'r');
    xlabel('Time relative to Mitosis (hrs)');
    ylabel('CDK2 activity'); ylabel('ERK activity');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]);
    saveas(gcf,'h:\Downloads\Fig2.jpg');
end
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heatmapoption=1;
if heatmapoption
    figure;
    heatmapmax=2.5;
    data2mother=saturatepertrace(data2mother);
    heatdatatotal=[data2mother*heatmapmax dataCdk2];
    heatdatatotal(heatdatatotal>heatmapmax)=heatmapmax;
    heatdatatotal(heatdatatotal<0)=0;
    heatdatatotal(isnan(heatdatatotal))=0; %make invalid data black
    traceid=1:numtraces;
    imagesc(time,traceid,heatdatatotal);
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
    saveas(gcf,'h:\Downloads\Fig3.jpg');
end
fprintf('%0.0f traces\n',numtraces);
end


function [tracesCdk2,traces2,samplestats]=main(datadir,shot,minlength)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(genealogy));
%samplecells=find(~isnan(genealogy) & ~isnan(tracedata(:,end,1))); %only cells that make it to end
%%% append mother data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
linkedtracedata=tracedata;
for i=1:numel(samplecells)
    prevframe=tracestats(samplecells(i),1)-1;
    linkedtracedata(samplecells(i),1:prevframe,:)=tracedata(genealogy(samplecells(i)),1:prevframe,:);
end
tracedata=linkedtracedata;
%%% gate traces based on signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals=tracedata(samplecells,:,5);
signals2=tracedata(samplecells,:,7)./tracedata(samplecells,:,5);
stats=tracestats(samplecells,:);
%hist(signals(:),0:10:2010); xlim([0 2000]);
minthresh=150; %HW:150
[samplecells,samplestats]=remove_short_neg_noisy_Cdk2(samplecells,signals,signals2,stats,minlength,minthresh);
%%% optional second gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signals=tracedata(samplecells,:,8);
%hist(signals(:),100);
minthresh=0;
[samplecells,samplestats]=remove_short_neg_noisy(samplecells,signals,samplestats,minlength,minthresh);
%%% calculate signal of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesCdk2=tracedata(:,:,7)./tracedata(:,:,5); %Sabrina:6/5 HeeWon:7/5
%%% optional second signal (just repeat tracesCdk2 if not wanted) %%%%%%%%%
traces2=tracedata(:,:,8)./tracedata(:,:,6);   %HeeWon:8/6
%%% smooth data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:numel(samplecells)
%     tracenan=isnan(tracesCdk2(samplecells(i),:));
%     tracesCdk2(samplecells(i),:)=smooth(tracesCdk2(samplecells(i),:));
%     traces2(samplecells(i),:)=smooth(traces2(samplecells(i),:));
%     tracesCdk2(samplecells(i),tracenan)=NaN;
%     traces2(samplecells(i),tracenan)=NaN;
% end
%%% remove non-sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesCdk2=tracesCdk2(samplecells,:);
traces2=traces2(samplecells,:);
end