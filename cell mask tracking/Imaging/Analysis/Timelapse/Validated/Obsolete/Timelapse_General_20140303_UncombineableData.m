function Timelapse_General(conditions,datadir)
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=1; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
if quiescentanalysis
    motheroption=1; daughteroption=2;
end
IFoption=0; %0:no IFdata 1:IFdata
drugspike=0;
minlengthtrace=15; %CDK2fatecall-minimum=15 Geminin-minimum=10
minlengthmother=10;
framesperhr=5;
ymin=0.3; ymax=2.5; ystring='CDK2 activity';
%ymin=0; ymax=1; ystring='Geminin';

condnum=size(conditions,1);
xstring='Time relative to Mitosis (hrs)';
%xstring='Time relative to drugspike (hrs)';
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    tracedata=[];
    tracestats=[];
    motherstats=[];
    IFdata=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [tracedatatemp,tracestatstemp,motherstatstemp,IFdatatemp]=gathertracedata_1(datadir,shot,motheroption,daughteroption,IFoption);
                tracedata=[tracedata;tracedatatemp];
                tracestats=[tracestats;tracestatstemp];
                motherstats=[motherstats;motherstatstemp];
                IFdata=[IFdata;IFdatatemp];
            end
        end
    end
end
%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucchannel=6; cytochannel=8; %6 8
maxthresh=125; %threshold above which max of each trace must be
noisethresh=0.5; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
[traces1,badtraces1]=gate_Cdk2_1(tracedata,nucchannel,cytochannel,tracestats,minlengthtrace,maxthresh,noisethresh,quiescentanalysis);

%%% gate pipdeg data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channelPIP=7;
% maxthresh=100;  %threshold above which max of each trace must be
% minthresh=10; %threshold below which min of each trace must be
% [traces1,badtraces1]=gate_lengthandrange(tracedata,tracestats,channelPIP,minlengthtrace,maxthresh,minthresh);

%%% gate Geminin data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channelGem=7;
% maxthresh=150;  %threshold above which max of each trace must be
% minthresh=100; %threshold below which min of each trace must be
% [traces2,badtraces2]=gate_lengthandrange(tracedata,tracestats,channelGem,minlengthtrace,maxthresh,minthresh);

%%% gate Erk data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channelErk=7;
% maxthresh=0;     %threshold above which max of each trace must be
% minthresh=10000; %threshold below which min of each trace must be
% [traces2,badtraces2]=gate_lengthandrange(tracedata,motherstats,channelErk,minlengthmother,maxthresh,minthresh);

%%% transform and gate IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IFoption
    channelIF=9;
    IFvals=IFdata(:,channelIF); badtracesIF=IFvals<1; IFvals=log2(IFvals);
else
    IFvals=ones(size(badtraces1))*NaN;
end
%%% combine gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%badtraces=zeros(size(tracedata,1),1);
badtraces=badtraces1;
%badtraces=badtraces1 & badtraces2;
%badtraces=badtraces1 | badtracesIF;
%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces1=traces1(~badtraces,:);
%traces2=traces2(~badtraces,:);
tracestats=tracestats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
IFvals=IFvals(~badtraces);
%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normoption=0;
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
%traces1=normalizetraces_3(traces1,tracestats,normoption);
%%% categorize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Cdk2inc,Cdk2low]=categorizeCdk2fate(traces1,tracestats,minlengthtrace);
%framereltomitosis=5; prctilethresh=[75 25]; [trace2high,trace2low]=categorizetracebytime(traces2,tracestats,framereltomitosis,prctilethresh);

%%% get Points Of Interest (absolute rather than relative to mitosis) %%%%%
[POI_Cdk2,badtraces_Cdk2]=GetPointsOfInterest_2('Cdk2',traces1,tracestats,minlengthtrace); %detect Cdk2start for daughtertraces
%[POI_Geminin,badtraces_Geminin]=GetPointsOfInterest_2('Geminin',traces2,tracestats,minlengthtrace);
%[POI_PipdegFall,badtraces_Pipdeg]=GetPointsOfInterest_2('PipdegFall',traces1,tracestats,minlengthtrace);
POIcalc=POI_Cdk2;
%POIcalc=tracestats(:,3);
%%% gate POI calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%badtracesallPOI=zeros(size(POIcalc,1),1); %no bad trace gating
badtracesallPOI=badtraces_Cdk2;
%badtracesallPOI=badtraces_Cdk2 | badtraces_Geminin;
POIcalc=POIcalc(~badtracesallPOI);
tracesgated=traces1(~badtracesallPOI,:); tracestats=tracestats(~badtracesallPOI,:); motherstats=motherstats(~badtracesallPOI,:);

%%% align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cdk2inc=Cdk2inc(~badtracesallPOI); Cdk2low=Cdk2low(~badtracesallPOI);
%alignPOI=tracestats(:,1); %align by mitosis
alignPOI=POIcalc;
%alignPOI=ones(size(tracesgated,1),1)*drugspike; %align by drugspike
[alignedtime,aligneddata]=aligntraces_3(tracesgated,alignPOI,tracestats,motherstats); alignedtime=alignedtime/framesperhr;
%%% plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesmoothoption=1; %1:smoothtraces 0:don't smooth
%plottraces_1(aligneddata,alignedtime,xstring,ystring,ylimits,tracesmoothoption); %plot traces without categorization
comparetrends_1(aligneddata,alignedtime,Cdk2inc,Cdk2low,xstring,ystring,[ymin ymax],tracesmoothoption);
%comparetrends_1(aligneddata,time,trace2high,trace2low,xstring,ystring,ylimits,tracesmoothoption);

%%% plot IF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if IFoption
    figure; hold on;
    IFtime=tracestats(:,3)/framesperhr;
    %scatter(IFtime(Cdk2inc),IFvals(Cdk2inc),'b','o'); scatter(IFtime(Cdk2low),IFvals(Cdk2low),'r','o','markerfacecolor','r');
    scatter(IFtime,IFvals,'b','o');
    xlabel(xstring); ylabel('mean(pRb(807/811))');
    set(gcf,'color','w','PaperPosition',[0 0 4 7]); %saveas(gcf,'h:\Downloads\FigIF.jpg');
end
%%% sort POIcalc for heatmap figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sortedPOIcalc,sortidx]=sort(POIcalc);
sortedaligneddata=aligneddata(sortidx,:);
%%% generate heatmap (POI must be relative to drugspike/mitosis/etc.) %%%%%
heatmapmin=ymin; heatmapmax=ymax;
%datasorted=saturatepertrace(datasorted);
markoption=1; %1:overlay Cdk2start 2:don't overlay marks
makeheatmaps_1(sortedaligneddata,sortedPOIcalc,alignedtime,heatmapmax,heatmapmin,markoption,xstring,framesperhr);
%%% record number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(aligneddata,1);
fprintf('%0.0f traces\n',numtraces);
end