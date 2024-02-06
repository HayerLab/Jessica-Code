function Align_Mitosis_CDK2(conditions,datadir)
alignoption=2; %1:mitosis 2:drugspike
signaloption=1; %1:Geminin 2:Signal2
motheroption=1;
minlength=40;
drugspike=40;
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltracesGeminin=[];
    alltracesSignal2=[];
    allstats=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [tracesGeminin,tracesSignal2,stats]=main(datadir,shot,minlength);
                alltracesGeminin=[alltracesGeminin;tracesGeminin];
                alltracesSignal2=[alltracesSignal2;tracesSignal2];
                allstats=[allstats;stats];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
[numtraces,numframes]=size(alltracesGeminin);
sampleid=(1:numtraces)';
detectrisetime=1;
if detectrisetime
    %%% detect onset of CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relative=1; %returns Cdk2start time relative to mitosis
    [GemRise,badtraces]=getGemininFeatures(alltracesGeminin,allstats,minlength,relative);
    GemRise(isnan(GemRise))=0;
    sampleid(badtraces>0)=[];
end
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sorttraces=1;
if sorttraces
    if alignoption==1
        [~,idx]=sort(GemRise(sampleid));
    elseif alignoption==2
        [~,idx]=sort(allstats(sampleid,1));
    end
    ordered=sampleid(idx);
    GemRise=GemRise(ordered);
    alltracesGeminin=alltracesGeminin(ordered,:);
    alltracesSignal2=alltracesSignal2(ordered,:);
    allstats=allstats(ordered,:);
end
%%% align data to mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltracesGeminin,1);
dataGeminin=ones(numtraces,numframes)*NaN;
dataSignal2=ones(numtraces,numframes)*NaN;
if alignoption==1
    for i=1:numtraces
        dataGeminin(i,1:allstats(i,3))=alltracesGeminin(i,allstats(i,1):allstats(i,2));
        dataSignal2(i,1:allstats(i,3))=alltracesSignal2(i,allstats(i,1):allstats(i,2));
    end
    maxlength=prctile(allstats(:,3),90);
    dataGeminin=dataGeminin(:,1:maxlength);
    dataSignal2=dataSignal2(:,1:maxlength);
elseif alignoption==2
    for i=1:numtraces
        dataGeminin(i,allstats(i,1):allstats(i,2))=alltracesGeminin(i,allstats(i,1):allstats(i,2));
        dataSignal2(i,allstats(i,1):allstats(i,2))=alltracesSignal2(i,allstats(i,1):allstats(i,2));
    end
    maxlength=prctile(allstats(:,2),90);
    dataGeminin=dataGeminin(:,1:maxlength);
    dataSignal2=dataSignal2(:,1:maxlength);
end
%%% align mother data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mothermax=0; dataGemininmother=[]; dataSignal2mother=[];
if motheroption && alignoption==1
    dataGemininmother=ones(numtraces,numframes)*NaN;
    dataSignal2mother=ones(numtraces,numframes)*NaN;
    for i=1:numtraces
        firstframe=find(~isnan(alltracesGeminin(i,:)),1,'first');
        lastframe=allstats(i,1)-1;
        duration=lastframe-firstframe+1;
        dataGemininmother(i,numframes-duration+1:numframes)=alltracesGeminin(i,firstframe:lastframe);
        dataSignal2mother(i,numframes-duration+1:numframes)=alltracesSignal2(i,firstframe:lastframe);
    end
    mothermax=80;
    dataGemininmother=dataGemininmother(:,numframes-mothermax+1:numframes);
    dataSignal2mother=dataSignal2mother(:,numframes-mothermax+1:numframes);
    %%% concatenate mother & daughter data %%%
    if signaloption==1
        datatotal=[dataGemininmother dataGeminin];
    elseif signaloption==2
        datatotal=[dataSignal2mother dataSignal2];
    end
elseif motheroption && alignoption==2
    dataGeminin=alltracesGeminin;
    dataSignal2=alltracesSignal2;
    if signaloption==1
        datatotal=dataGeminin;
    elseif signaloption==2
        datatotal=dataSignal2;
    end
end
if alignoption==1
    time=(-mothermax+1:maxlength)/framesperhr;
elseif alignoption==2
    time=(1-drugspike:maxlength-drugspike)/framesperhr;
end
%%% generate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traceoption=0;
if traceoption
    figure, hold on;
    for i=1:size(datatotal,1)
        line(time,datatotal(i,:),'color','b');
    end
    xlabel('Time relative to Mitosis (hrs)');
    ylabel('Geminin activity'); %ylabel('Signal2 activity');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]);
    saveas(gcf,'h:\Downloads\Fig1.jpg');
    
    figure, hold on;
    binstep=1;
    bincurveshade(time,dataCdk2inc,binstep,'b');
    bincurveshade(time,dataCdk2low,binstep,'r');
    xlabel('Time relative to Mitosis (hrs)');
    ylabel('Geminin activity'); ylabel('Signal2 activity');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]);
    saveas(gcf,'h:\Downloads\Fig2.jpg');
end
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heatmapoption=1;
if heatmapoption
    figure;
    heatdatatotal=datatotal;
    heatmapmax=1; heatmapmin=0;
    heatdatatotal(heatdatatotal>heatmapmax)=heatmapmax;
    heatdatatotal(heatdatatotal<heatmapmin)=heatmapmin;
    heatdatatotal(isnan(heatdatatotal))=heatmapmin; %make invalid data black
    traceid=1:numtraces;
    imagesc(time,traceid,heatdatatotal);
    cmap=colormap(jet);
    cmap(1,:)=[0 0 0];
    colormap(cmap);
    %%% overlay markers at timepoints of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    addmarks=1;
    if addmarks && alignoption==1
        hold on;
        onlyreal=find(GemRise~=0);
        GemRise_onlyreal=GemRise(onlyreal);
        traceid_onlyreal=traceid(onlyreal);
        scatter((GemRise_onlyreal)/framesperhr,traceid_onlyreal,8,'r','*');
    elseif addmarks && alignoption==2
        hold on;
        scatter(zeros(numtraces,1),1:numtraces,12,'r','*');
        %scatter((allstats(:,1)-drugspike)/framesperhr,1:numtraces,12,'g','*');
    end
    %%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xlabel('Time relative to Mitosis (hrs)'); ylabel('Individual traces');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]);
    saveas(gcf,'h:\Downloads\Fig3.jpg');
end
fprintf('%0.0f traces\n',numtraces);
end


function [Geminin,Signal2,daughterstats]=main(datadir,shot,minlength)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(genealogy));
%%% get full durations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fulltracedata=linkancestry(tracedata,tracestats,samplecells);
%%% append mother data and record track mother stats %%%%%%%%%%%%%%%%%%%%%%
motherstats=ones(size(tracestats,1),size(tracestats,2)+1)*NaN;
for i=1:numel(samplecells)
    s=samplecells(i);
    m=tracestats(s,4);
    ancestrylength=tracestats(m,2)-find(~isnan(fulltracedata(s,:,1)),1,'first')+1;
    motherstats(s,:)=[tracestats(m,:) ancestrylength];
    motherframes=motherstats(s,1):motherstats(s,2);
    tracedata(s,motherframes,:)=tracedata(m,motherframes,:);
end
tracedata=tracedata(samplecells,:,:);
daughterstats=tracestats(samplecells,:);
motherstats=motherstats(samplecells,:);
%%% smooth data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
Geminin=tracedata(:,:,6);
Signal2=tracedata(:,:,5);
for i=1:numtraces
    realframes=find(~isnan(Geminin(i,:)));
    Geminin(i,realframes)=smooth(Geminin(i,realframes));
    Signal2(i,realframes)=smooth(Signal2(i,realframes));
end
%%% gate traces based on Geminin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(min(Geminin,[],2),0:10:1010); xlim([0 1000]);
%hist(max(Geminin,[],2),0:10:1010); xlim([0 1000]);
maxthresh=100;
minthresh=1000;
badtracesGeminin=gate_highlow(Geminin,daughterstats,minlength,maxthresh,minthresh);
%%% gate traces based on pipdeg sensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(max(tracespipdeg,[],2),0:20:2020); xlim([0 2000]);
%hist(min(tracespipdeg,[],2),0:2:202); xlim([0 200]);
maxthresh=100; %50prc:50 10nuc50prc:100-->200
minthresh=1000; %50prc:25 10nuc50prc:100-->75
badtracesSignal2=gate_highlow(Signal2,daughterstats,minlength,maxthresh,minthresh); %gate length by mother or daughter recs
badtraces=badtracesGeminin | badtracesSignal2;
%%% remove non-sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Geminin=Geminin(~badtraces,:);
Signal2=Signal2(~badtraces,:);
daughterstats=daughterstats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
%%% normalize pipdeg traces by max in mother %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(Geminin,1)
    maxvalGeminin=max(Geminin(i,motherstats(i,1):motherstats(i,2)));
    Geminin(i,:)=Geminin(i,:)/maxvalGeminin;
    maxvalSignal2=max(Signal2(i,motherstats(i,1):motherstats(i,2)));
    Signal2(i,:)=Signal2(i,:)/maxvalSignal2;    
end
end