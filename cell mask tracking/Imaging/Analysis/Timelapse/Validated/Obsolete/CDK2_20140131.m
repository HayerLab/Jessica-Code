function CDK2(conditions,datadir)
alignoption=2; %1:mitosis 2:drugspike
signaloption=1; %1:CDK2 2:Second Signal
drugspike=10;
minlengthdaughter=20; %minimum=15
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
    tracedata=[];
    daughterstats=[];
    motherstats=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                %shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [tracedatatemp,daughterstatstemp,motherstatstemp]=getdatadaughtermother(datadir,shot);
                tracedata=[tracedata;tracedatatemp];
                daughterstats=[daughterstats;daughterstatstemp];
                motherstats=[motherstats;motherstatstemp];
            end
        end
    end
end
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
%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesCdk2=tracesCdk2(~badtraces,:);
traces2=traces2(~badtraces,:);
daughterstats=daughterstats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
%%% normalize traces by max in daughter or mother %%%%%%%%%%%%%%%%%%%%%%%%%
%traces2=normalizetraces(traces2,motherstats);
%%% categorize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Cdk2inc,Cdk2low]=categorizeCdk2fate(tracesCdk2,daughterstats,minlengthdaughter);
%framereltomitosis=5; prctilethresh=[75 25]; [trace2high,trace2low]=categorizetracebytime(traces2,daughterstats,framereltomitosis,prctilethresh);
%%% align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if signaloption==1
    alltraces=tracesCdk2;
    ystring='CDK2 activity'; ylimits=[0 2.5];
elseif signaloption==2
    alltraces=traces2;
    ystring='ERK activity'; ylimits=[0 5];
end
historyoption=1; %1:motherdata 2:allancestry
motherperc=80; %mother length to display (percentile of durations)
daughterperc=80; %daughter length to display
[time,datatotal]=alignmotherdaughter2(alltraces,daughterstats,motherstats,alignoption,drugspike,historyoption,motherperc,daughterperc,framesperhr);
%%% plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottraces(datatotal,time,xstring,ystring,ylimits); %plot traces without categorization
comparetrends(datatotal,time,Cdk2inc,Cdk2low,xstring,ystring,ylimits);
%comparetrends(datatotal,time,trace2high,trace2low,xstring,ystring,ylimits);
%%% get and sort Cdk2start times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[datasorted,Cdk2start]=sortCdk2starts(tracesCdk2,datatotal,daughterstats,minlengthdaughter,alignoption); %detect Cdk2start for daughtertraces
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
numtraces=size(datatotal,1);
fprintf('%0.0f traces\n',numtraces);
end