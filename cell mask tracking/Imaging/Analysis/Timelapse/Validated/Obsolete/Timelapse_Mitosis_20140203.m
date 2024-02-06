function Timelapse_Mitosis(conditions,datadir)
IFoption=1;
alignoption=1; %1:mitosis 2:drugspike
drugspike=10;
minlengthdaughter=15; %CDK2fatecall-minimum=15
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
    IFdata=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                %shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [tracedatatemp,daughterstatstemp,motherstatstemp,IFdatatemp]=getdatadaughtermother_test(datadir,shot,IFoption);
                tracedata=[tracedata;tracedatatemp];
                daughterstats=[daughterstats;daughterstatstemp];
                motherstats=[motherstats;motherstatstemp];
                IFdata=[IFdata;IFdatatemp];
            end
        end
    end
end
%%% smooth and gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucchannel=5; cytochannel=7;
maxthresh=150; %threshold above which max of each trace must be
noisethresh=0.15; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
[traces1,badtraces1]=gate_Cdk2(tracedata,nucchannel,cytochannel,daughterstats,minlengthdaughter,maxthresh,noisethresh);

%%% smooth and gate Geminin data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channel1=6;
maxthresh=150;  %threshold above which max of each trace must be
minthresh=100; %threshold below which min of each trace must be
[traces2,badtraces2]=gate_smoothlengthrange(tracedata,daughterstats,channel1,minlengthdaughter,maxthresh,minthresh);

%%% gate Erk data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% traces2gating=tracedata(:,:,5);
% traces2=tracedata(:,:,7)./tracedata(:,:,5);
% maxthresh=0;     %threshold above which max of each trace must be
% minthresh=10000; %threshold below which min of each trace must be
% badtraces2=gate_lengthrange(traces2gating,motherstats,minlengthmother,maxthresh,minthresh);
% [traces2,badtraces2]=gate_smoothlengthrange(tracedata,daughterstats,channel7,minlengthmother,maxthresh,minthresh);

%%% transform and gate IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelIF=8;
IFvals=IFdata(:,channelIF); IFvals(IFvals<1)=1; IFvals=log2(IFvals);

%%% combine gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%badtraces=badtraces1;
badtraces=badtraces1 | badtraces2;
%badtraces=badtraces1;
%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces1=traces1(~badtraces,:);
traces2=traces2(~badtraces,:);
daughterstats=daughterstats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
IFvals=IFvals(~badtraces);
%%% normalize traces by max in daughter or mother %%%%%%%%%%%%%%%%%%%%%%%%%
%traces1=normalizetraces(traces1,motherstats);
%%% categorize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Cdk2inc,Cdk2low]=categorizeCdk2fate(traces1,daughterstats,minlengthdaughter);
%framereltomitosis=5; prctilethresh=[75 25]; [trace2high,trace2low]=categorizetracebytime(traces2,daughterstats,framereltomitosis,prctilethresh);
%%% align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
historyoption=1; %1:motherdata 2:allancestry
motherperc=90; %mother length to display (percentile of durations)
daughterperc=90; %daughter length to display
[alignedtime,aligneddata]=alignmotherdaughter2(traces1,daughterstats,motherstats,alignoption,drugspike,historyoption,motherperc,daughterperc,framesperhr);
%%% plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin=0.3; ymax=2.5; ystring='CDK2 activity';
%ymin=0; ymax=1; ystring='Geminin';
ylimits=[ymin ymax];
%plottraces(aligneddata,time,xstring,ystring,ylimits); %plot traces without categorization
comparetrends(aligneddata,alignedtime,Cdk2inc,Cdk2low,xstring,ystring,ylimits);
%comparetrends(aligneddata,time,trace2high,trace2low,xstring,ystring,ylimits);
%%% plot IF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFtime=daughterstats(:,3)/framesperhr;
figure; hold on;
scatter(IFtime(Cdk2inc),IFvals(Cdk2inc),'b','o');
scatter(IFtime(Cdk2low),IFvals(Cdk2low),'r','o','markerfacecolor','r');
xlabel(xstring); ylabel('IF levels'); ylim([7 15]); % 7 15
%%% get and sort Cdk2start times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[datasorted,POI]=GetPointsOfInterest('Cdk2',traces1,aligneddata,daughterstats,minlengthdaughter,alignoption); %detect Cdk2start for daughtertraces
%[datasorted,POI]=GetPointsOfInterest('Geminin',traces1,aligneddata,daughterstats,minlengthdaughter,alignoption);
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heatmapmin=ymin; heatmapmax=ymax;
%datasorted=saturatepertrace(datasorted);
markoption=1; %1:overlay Cdk2start 2:don't overlay marks
makeCdk2heatmaps(datasorted,POI,alignedtime,heatmapmax,heatmapmin,markoption,alignoption,xstring,framesperhr);
%%% record number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(aligneddata,1);
fprintf('%0.0f traces\n',numtraces);
end