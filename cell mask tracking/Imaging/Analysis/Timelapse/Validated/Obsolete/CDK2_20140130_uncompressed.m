function CDK2(conditions,datadir)
alignoption=1; %1:mitosis 2:drugspike
signaloption=1; %1:CDK2 2:Second Signal
drugspike=0;
minlengthdaughter=30;
minlengthmother=10;
framesperhr=5;
condnum=size(conditions,1);
if alignoption==1
    xstring='Time relative to Mitosis (hrs)';
elseif alignoption==2
    xstring='Time relative to drugspike (hrs)';
end
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltracesCdk2=[];
    alltraces2=[];
    alldaughterstats=[];
    allmotherstats=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [tracesCdk2,tracespipdeg,daughterstats,motherstats]=main(datadir,shot,minlengthdaughter,minlengthmother);
                alltracesCdk2=[alltracesCdk2;tracesCdk2];
                alltraces2=[alltraces2;tracespipdeg];
                alldaughterstats=[alldaughterstats;daughterstats];
                allmotherstats=[allmotherstats;motherstats];
            end
        end
    end
end
%%% categorize daughter CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltracesCdk2,1);
earlyval=ones(numtraces,1)*NaN;
maxval=ones(numtraces,1)*NaN;
%lastval=ones(numtraces,1)*NaN;
for i=1:numtraces
    earlyval(i)=alltracesCdk2(i,alldaughterstats(i,1)+minlengthdaughter-1); %Sabrina uses 2hr
    maxval(i)=max(alltracesCdk2(i,alldaughterstats(i,1)+minlengthdaughter-1:alldaughterstats(i,2))); %max from minlength-end
    %lastval(i)=alltracesCdk2(i,alldaughterstats(i,2));
end
%figure,hist(earlyval,0.3:0.02:3);
mincutoff=0.60;
maxcutoff=0.80;
Cdk2inc=earlyval>mincutoff & maxval>maxcutoff+0.2;
Cdk2low=earlyval<mincutoff & maxval<maxcutoff-0.2;
%%% categorize alltraces2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val2=ones(numtraces,1)*NaN;
postmitoticframe=5;
for i=1:numtraces
    val2(i)=alltraces2(i,alldaughterstats(i,1)+postmitoticframe-1);
end
%figure,hist(val2,100);
highval2=val2>prctile(val2,75);
lowval2=val2<prctile(val2,25);
%%% align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if signaloption==1
    alltraces=alltracesCdk2;
    ystring='CDK2 activity'; ylimits=[0 2.5];
elseif signaloption==2
    alltraces=alltraces2;
    ystring='ERK activity'; ylimits=[0 5];
end
historyoption=1; %1:motherdata 2:allancestry
motherperc=80;
daughterperc=80;
[time,datatotal]=alignmotherdaughter2(alltraces,alldaughterstats,allmotherstats,alignoption,historyoption,motherperc,daughterperc,framesperhr);
%%% plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%comparetrends(datatotal,time,Cdk2inc,Cdk2low,xstring,ystring,ylimits);
comparetrends(datatotal,time,highval2,lowval2,xstring,ystring,ylimits);

%%% get and sort Cdk2start times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[datasorted,Cdk2start]=sortCdk2starts(alltracesCdk2,datatotal,alldaughterstats,minlengthdaughter,alignoption);
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if signaloption==1
    heatmapmax=2.5; heatmapmin=0.3;
elseif signaloption==2
    datasorted=saturatepertrace(datasorted);
    heatmapmax=1; heatmapmin=0;
end
markoption=1; %1:overlay Cdk2start 2:don't overlay marks
makeCdk2heatmaps(datasorted,Cdk2start,time,heatmapmax,heatmapmin,markoption,alignoption,xstring,framesperhr);
%%% record number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%0.0f traces\n',numtraces);
end

function [tracesCdk2,traces2,daughterstats,motherstats]=main(datadir,shot,minlengthdaughter,minlengthmother)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(genealogy));
%%% get full durations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracedata=linkancestry(tracedata,tracestats,samplecells);
%%% record mother stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motherstats=ones(size(tracestats,1),size(tracestats,2)+1)*NaN;
for i=1:numel(samplecells)
    s=samplecells(i);
    m=tracestats(s,4);
    ancestrylength=tracestats(m,2)-find(~isnan(tracedata(s,:,1)),1,'first')+1;
    motherstats(s,:)=[tracestats(m,:) ancestrylength];
end
tracedata=tracedata(samplecells,:,:);
daughterstats=tracestats(samplecells,:);
motherstats=motherstats(samplecells,:);
%%% smooth CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
DHBnuc=tracedata(:,:,8);
tracesCdk2=tracedata(:,:,8)./tracedata(:,:,6);
for i=1:numtraces
    realframes=find(~isnan(DHBnuc(i,:)));
    DHBnuc(i,realframes)=smooth(DHBnuc(i,realframes));
    tracesCdk2(i,realframes)=smooth(tracesCdk2(i,realframes));
end
%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(DHBnuc(:),0:10:2010); xlim([0 2000]);
%hist(max(DHBnuc,[],2),0:20:1020); xlim([0 1000]);
thresh=150; %HW:150 pipdeg0715/0719:localbg=75 sigbg=75
%badtracesCdk2=gate_Cdk2(DHBnuc,tracesCdk2,daughterstats,minlengthdaughter,thresh);
noisethresh=0.5; %10x:0.15  20x:0.5;
badtracesCdk2=gate_Cdk2_noisecalcdiff1(DHBnuc,tracesCdk2,daughterstats,minlengthdaughter,thresh,noisethresh);
%%% gate traces based on pipdeg sensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces2gating=tracedata(:,:,7);
%traces2=tracedata(:,:,6);
traces2=tracedata(:,:,7)./tracedata(:,:,5);
%hist(min(traces2gating,[],2),0:20:2020); xlim([0 2000]);
%hist(min(traces2gating,[],2),0:80:8080); xlim([0 8000]);
%hist(log2(min(traces2gating,[],2)),0:0.1:15.1); xlim([0 15]);
maxthresh=0;     %HeeWon:0
minthresh=10000; %HeeWon:10000
badtraces2=gate_lengthrange(traces2gating,motherstats,minlengthmother,maxthresh,minthresh);
badtraces=badtracesCdk2 | badtraces2;
%%% remove non-sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesCdk2=tracesCdk2(~badtraces,:);
traces2=traces2(~badtraces,:);
daughterstats=daughterstats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
%%% normalize pipdeg traces by max in mother %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalizeoption=0;
if normalizeoption
    for i=1:size(traces2,1)
        maxval=max(traces2(i,motherstats(i,1):motherstats(i,2)));
        traces2(i,:)=traces2(i,:)/maxval;
    end
end
end