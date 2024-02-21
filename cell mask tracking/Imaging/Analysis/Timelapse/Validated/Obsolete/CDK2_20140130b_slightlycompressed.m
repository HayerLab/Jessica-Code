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
                shot=wellnum2str(row,col,site);
                %shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
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
comparetrends(datatotal,time,Cdk2inc,Cdk2low,xstring,ystring,ylimits);
%comparetrends(datatotal,time,highval2,lowval2,xstring,ystring,ylimits);
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
[daughterstats,motherstats,tracedata]=getmotherstats(tracedata,tracestats,samplecells);

%%% smooth and gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucchannel=5; cytochannel=7;
maxthresh=150; %threshold above which max of each trace must be
noisethresh=0.15; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
[tracesCdk2,badtracesCdk2]=gate_Cdk2(tracedata,nucchannel,cytochannel,daughterstats,minlengthdaughter,maxthresh,noisethresh);

%%% gate traces based on second sensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces2gating=tracedata(:,:,5);
traces2=tracedata(:,:,7)./tracedata(:,:,5);
maxthresh=0;     %threshold above which max of each trace must be
minthresh=10000; %threshold below which min of each trace must be
badtraces2=gate_lengthrange(traces2gating,motherstats,minlengthmother,maxthresh,minthresh);

%%% combine gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces=badtracesCdk2 | badtraces2;
%%% remove non-sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesCdk2=tracesCdk2(~badtraces,:);
traces2=traces2(~badtraces,:);
daughterstats=daughterstats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
%%% normalize traces by max in daughter or mother %%%%%%%%%%%%%%%%%%%%%%%%%%
traces2=normalizetraces(traces2,motherstats);
end