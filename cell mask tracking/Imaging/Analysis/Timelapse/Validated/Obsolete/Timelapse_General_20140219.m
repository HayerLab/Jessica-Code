function Timelapse_Mitosis(conditions,datadir)
motheroption=1; %0:no gating 1:mothers 2:no mothers
daughteroption=1: %0:no gating 1:daughters 2:no daughters
IFoption=0; %0:no IFdata 1:IFdata
alignoption=1; %1:mitosis 2:drugspike
drugspike=0;
minlengthdaughter=15; %CDK2fatecall-minimum=15 Geminin-minimum=10
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
                [tracedatatemp,daughterstatstemp,motherstatstemp,IFdatatemp]=gathertracedata1(datadir,shot,motheroption,daughteroption,IFoption);
                tracedata=[tracedata;tracedatatemp];
                daughterstats=[daughterstats;daughterstatstemp];
                motherstats=[motherstats;motherstatstemp];
                IFdata=[IFdata;IFdatatemp];
            end
        end
    end
end
%%% smooth and gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucchannel=6; cytochannel=8;
maxthresh=150; %threshold above which max of each trace must be
noisethresh=0.15; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
[traces1,badtraces1]=gate_Cdk2(tracedata,nucchannel,cytochannel,daughterstats,minlengthdaughter,maxthresh,noisethresh);

%%% smooth and gate Geminin data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channel1=7;
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
if IFoption
    channelIF=8;
    IFvals=IFdata(:,channelIF); badtracesIF=IFvals<1; IFvals=log2(IFvals);
else
    IFvals=ones(size(badtraces1))*NaN;
end
%%% combine gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%badtraces=badtraces1;
badtraces=badtraces1 & badtraces2;
%badtraces=badtraces1 | badtraces2 | badtracesIF;
%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces1=traces1(~badtraces,:);
traces2=traces2(~badtraces,:);
daughterstats=daughterstats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
IFvals=IFvals(~badtraces);
%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces2=normalizetraces(traces2);
%%% categorize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Cdk2inc,Cdk2low]=categorizeCdk2fate(traces1,daughterstats,minlengthdaughter);
%framereltomitosis=5; prctilethresh=[75 25]; [trace2high,trace2low]=categorizetracebytime(traces2,daughterstats,framereltomitosis,prctilethresh);
%%% align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motherperc=90; %mother length to display (percentile of durations)
daughterperc=90; %daughter length to display
[alignedtime,aligneddata]=aligntraces(traces1,daughterstats,motherstats,alignoption,drugspike,mitosisoption,motherperc,daughterperc,framesperhr);
%%% plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin=0.3; ymax=2.5; ystring='CDK2 activity';
%ymin=0; ymax=1; ystring='Geminin';
ylimits=[ymin ymax];
%plottraces(aligneddata,alignedtime,xstring,ystring,ylimits); %plot traces without categorization
comparetrends(aligneddata,alignedtime,Cdk2inc,Cdk2low,xstring,ystring,ylimits);
%comparetrends(aligneddata,time,trace2high,trace2low,xstring,ystring,ylimits);
%%% get and sort Points Of Interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sortedidx,sortedPOI]=GetPointsOfInterest('Cdk2',traces1,daughterstats,minlengthdaughter,alignoption); %detect Cdk2start for daughtertraces
%[sortedidx,sortedPOI]=GetPointsOfInterest('Geminin',traces2,daughterstats,minlengthdaughter,alignoption);
sortedaligneddata=aligneddata(sortedidx,:);
%%% plot IF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IFoption
    figure; hold on;
    %IFtime=daughterstats(:,3)/framesperhr; scatter(IFtime(Cdk2inc),IFvals(Cdk2inc),'b','o'); scatter(IFtime(Cdk2low),IFvals(Cdk2low),'r','o','markerfacecolor','r');
    IFtime=(daughterstats(sortedidx,3)-sortedPOI)/framesperhr; IFtime(isnan(IFtime))=0; scatter(IFtime,IFvals(sortedidx),'b','o');
    xlabel(xstring); ylabel('CycA Levels');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]); saveas(gcf,'h:\Downloads\FigIF.jpg');
end
%%% generate heatmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heatmapmin=ymin; heatmapmax=ymax;
%datasorted=saturatepertrace(datasorted);
markoption=1; %1:overlay Cdk2start 2:don't overlay marks
makeheatmaps(sortedaligneddata,sortedPOI,alignedtime,heatmapmax,heatmapmin,markoption,alignoption,xstring,framesperhr);
%%% record number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(aligneddata,1);
fprintf('%0.0f traces\n',numtraces);
end