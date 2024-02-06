function Geminin(conditions,datadir)
signaloption=1; %1:CDK2 2:pipdeg
minlengthdaughter=30;
minlengthmother=10;
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltracesGeminin=[];
    alltracesSignal2=[];
    alldaughterstats=[];
    allmotherstats=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [tracesGeminin,traceSignal2,daughterstats,motherstats]=main(datadir,shot,minlengthdaughter,minlengthmother);
                alltracesGeminin=[alltracesGeminin;tracesGeminin];
                alltracesSignal2=[alltracesSignal2;traceSignal2];
                alldaughterstats=[alldaughterstats;daughterstats];
                allmotherstats=[allmotherstats;motherstats];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
[numtraces,numframes]=size(alltracesGeminin);
%%% extract pipdeg features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relative=0;
[degstarts,degends]=getdegstartandend_mother(alltracesSignal2,allmotherstats,minlengthmother,relative);
%%% sort traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sorttraces=1;
if sorttraces
    if alignoption==1
        [dummy,idx]=sort(Cdk2start(sampleid));
    elseif alignoption==2
        [dummy,idx]=sort(allstats(sampleid,1));
    end
    ordered=sampleid(idx);
    Cdk2start=Cdk2start(ordered);
    alltracesCdk2=alltracesCdk2(ordered,:);
    alltraces2=alltraces2(ordered,:);
    allstats=allstats(ordered,:);
end
%%% align daughter data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltracesGeminin,1);
dataGeminindaughter=ones(numtraces,numframes)*NaN;
dataSignal2daughter=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    dataGeminindaughter(i,1:alldaughterstats(i,3))=alltracesGeminin(i,alldaughterstats(i,1):alldaughterstats(i,2));
    dataSignal2daughter(i,1:alldaughterstats(i,3))=alltracesSignal2(i,alldaughterstats(i,1):alldaughterstats(i,2));
end
maxlength=prctile(alldaughterstats(:,3),90);
dataGeminindaughter=dataGeminindaughter(:,1:maxlength);
dataSignal2daughter=dataSignal2daughter(:,1:maxlength);
%%% align mother data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataGemininmother=ones(numtraces,numframes)*NaN;
dataSignal2mother=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    dataGemininmother(i,numframes-allmotherstats(i,3)+1:numframes)=alltracesGeminin(i,allmotherstats(i,1):allmotherstats(i,2));
    dataSignal2mother(i,numframes-allmotherstats(i,3)+1:numframes)=alltracesSignal2(i,allmotherstats(i,1):allmotherstats(i,2));
end
mothermax=80;
dataGemininmother=dataGemininmother(:,numframes-mothermax+1:numframes);
dataSignal2mother=dataSignal2mother(:,numframes-mothermax+1:numframes);
%%% concatenate mother & daughter data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if signaloption==1
    datatotal=[dataGemininmother dataGeminindaughter];
elseif signaloption==2
    datatotal=[dataSignal2mother dataSignal2daughter];
end
time=(-mothermax+1:maxlength)/framesperhr;
%%% generate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traceoption=1;
if traceoption
    figure, hold on;
    for i=1:size(datatotal,1)
        line(time,datatotal(i,:),'color','b');
    end
    xlabel('Time relative to Mitosis (hrs)');
    ylabel('Geminin');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]);
    saveas(gcf,'h:\Downloads\Fig1.jpg');
    
    figure, hold on;
    binstep=1;
    bincurveshade(time,datatotal,binstep,'b');
    xlabel('Time relative to Mitosis (hrs)');
    ylabel('Geminin'); ylabel('pipdeg activity');
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
        onlyreal=find(Cdk2start~=0);
        Cdk2start_onlyreal=Cdk2start(onlyreal);
        traceid_onlyreal=traceid(onlyreal);
        scatter((Cdk2start_onlyreal)/framesperhr,traceid_onlyreal,8,'r','*');
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


function [Geminin,Signal2,daughterstats,motherstats]=main(datadir,shot,minlengthdaughter,minlengthmother)
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
Geminin=tracedata(:,:,5);
Signal2=tracedata(:,:,6);
for i=1:numtraces
    realframes=find(~isnan(Geminin(i,:)));
    Geminin(i,realframes)=smooth(Geminin(i,realframes));
    Signal2(i,realframes)=smooth(Signal2(i,realframes));
end
%%% gate traces based on DHB sensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(DHBnuc(:),0:10:2010); xlim([0 2000]);
%hist(min(DHBnuc,[],2),0:20:1020); xlim([0 1000]);
maxthresh=75;
minthresh=50;
badtracesGeminin=gate_highlow(Geminin,daughterstats,minlengthdaughter,maxthresh,minthresh);
%%% gate traces based on pipdeg sensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(max(tracespipdeg,[],2),0:20:2020); xlim([0 2000]);
%hist(min(tracespipdeg,[],2),0:2:202); xlim([0 200]);
maxthresh=150; %50prc:50 10nuc50prc:100-->200
minthresh=75; %50prc:25 10nuc50prc:100-->75
badtracesSignal2=gate_highlow(Signal2,daughterstats,minlengthdaughter,maxthresh,minthresh); %gate length by mother or daughter recs
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