function Timelapse_General(conditions,datadir)
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=1; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
if quiescentanalysis
    motheroption=2; daughteroption=2;
end
IFoption=1; %0:no IFdata 1:IFdata
drugspike=1; %3C_CDK4iWee1i:14 %CyclingCDK4iWee1iCdc25i:59
minlengthtrace=11; %CDK2fatecall-minimum=11 Geminin-minimum=10 Pipdeg:5
minlengthmother=0; %PipDeg:50
framesperhr=5;
ymin=0.3; ymax=2.5; ystring='CDK2 activity';
%ymin=0; ymax=1; ystring='Geminin';
%ymin=0; ymax=1; ystring='p21dCy1';

%condnum=size(conditions,1);
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);

%xstring='Time relative to Mitosis (hrs)';
%xstring='Time relative to drugspike (hrs)';
xstring='Time since serum release (hrs)';
for i=1:uniquecondnum
    condrow=find(ismember(conditions(:,1),uniquenames{i}));
    tracedata=[];
    tracestats=[];
    motherstats=[];
    IFdata=[];
    for c=condrow'
        rowmat=cell2mat(conditions(c,2));
        colmat=cell2mat(conditions(c,3));
        sitemat=cell2mat(conditions(c,4));
        for row=rowmat
            for col=colmat
                for site=sitemat
                    %shot=wellnum2str(row,col,site);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    %shot=[num2str(row),'_',num2str(col),'_',num2str(site),'_100'];
                    [tracedatatemp,tracestatstemp,motherstatstemp,IFdatatemp]=gathertracedata_1(datadir,shot,motheroption,daughteroption,IFoption);
                    tracedata=[tracedata;tracedatatemp];
                    tracestats=[tracestats;tracestatstemp];
                    motherstats=[motherstats;motherstatstemp];
                    IFdata=[IFdata;IFdatatemp];
                end
            end
        end
    end
end

%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucchannel=6; cytochannel=7; %6 8
maxthresh=200; %threshold above which max of each trace must be
noisethresh=0.2; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
[traces1,badtraces1]=gate_Cdk2_1(tracedata,nucchannel,cytochannel,tracestats,minlengthtrace,maxthresh,noisethresh,quiescentanalysis);

%%% gate pipdeg data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channelPIP=7; %7
% maxthresh=50;  %threshold above which max of each trace must be %50
% maxpos=1; %0:anywhere 1:firstframe 2:lastframe
% minthresh=5; %threshold below which min of each trace must be %50
% minpos=1; %0:anywhere 1:mothertrace 2:daughtertrace
% %[traces1,badtraces1]=gate_lengthandrange(tracedata,tracestats,channelPIP,minlengthtrace,maxthresh,minthresh);
% %[traces1,badtraces1]=gate_lengthandrange_maxpos(tracedata,tracestats,channelPIP,minlengthtrace,maxthresh,minthresh,maxpos);
% [traces1,badtraces1]=gate_lengthandrange_minmother(tracedata,tracestats,motherstats,channelPIP,minlengthtrace,minlengthmother,maxthresh,minthresh,maxpos,minpos);

%%% gate H2B data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channelH2B=5;
% traces2=tracedata(:,:,channelH2B);

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
IFvals=ones(size(badtraces1))*NaN;
if IFoption
    channelIF=9;
    IFvals=IFdata(:,channelIF); badtracesIF=IFvals<1; IFvals(IFvals<1)=1; IFvals=log2(IFvals);
    Hoechst=IFdata(:,8).*IFdata(:,3); %Hoechst(Hoechst<1)=1; Hoechst=log2(Hoechst);
    EdU=IFdata(:,10).*IFdata(:,3); EdU(EdU<1)=1; EdU=log2(EdU);
    G1cells=Hoechst>2500000 & Hoechst<4500000 & EdU>10 & EdU<16;
    Scells=Hoechst>2500000 & Hoechst<7000000 & EdU>18 & EdU<24;
    %badtracesIF=IFvals1<=0 | IFvals2<=0;
    badtracesIF=(~G1cells) | badtracesIF;
end
%%% combine gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%badtraces=zeros(size(tracedata,1),1);
%badtraces=badtraces1;
%badtraces=badtraces1 & badtraces2;
badtraces=badtraces1 | badtracesIF;
%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces1=traces1(~badtraces,:);
%traces2=traces2(~badtraces,:);
tracestats=tracestats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
IFvals=IFvals(~badtraces);
%traces1=traces1./traces2;
%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normoption=1;
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
%traces2=normalizetraces_3(traces2,tracestats,normoption);
%%% categorize traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Cdk2inc,Cdk2low]=categorizeCdk2fate_2(traces1,tracestats,minlengthtrace);
%framereltomitosis=5; prctilethresh=[75 25]; [trace2high,trace2low]=categorizetracebytime(traces2,tracestats,framereltomitosis,prctilethresh);
%%% get Points Of Interest (absolute rather than relative to mitosis) %%%%%
%[POI_Cdk2,badtraces_Cdk2]=GetPointsOfInterest_2('Cdk2',traces1,tracestats,minlengthtrace); %detect Cdk2start for daughtertraces
%[POI_Geminin,badtraces_Geminin]=GetPointsOfInterest_2('Geminin',traces2,tracestats,minlengthtrace);
%[POI_PipdegFall,badtraces_Pipdeg]=GetPointsOfInterest_2('PipdegFall',traces1,tracestats,minlengthtrace);
%POIcalc=POI_Cdk2;
%POIcalc=POI_Cdk2-tracestats(:,1)+1;
%POIcalc=POI_Cdk2+timeoffset-drugspike;
POIcalc=tracestats(:,1)-drugspike;
%POIcalc=POI_PipdegFall-tracestats(:,1);
%%% determine POI gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtracesallPOI=zeros(size(POIcalc,1),1); %no bad trace gating
%badtracesallPOI=badtraces_Cdk2;
%badtracesallPOI=badtraces_Cdk2 | POIcalc>-15;
%badtracesallPOI=badtraces_Cdk2 | POIcalc>drugspike-1 | POIcalc<drugspike-20 | isnan(POIcalc);
%badtracesallPOI=badtraces_Cdk2 | POIcalc>timeoffset+5 | isnan(POIcalc);
%badtracesallPOI=badtraces_Cdk2 | badtraces_Geminin;
%badtracesallPOI=badtraces_Pipdeg;
%%% gate out badPOI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
POIcalc=POIcalc(~badtracesallPOI);
tracesgated=traces1(~badtracesallPOI,:); tracestats=tracestats(~badtracesallPOI,:); motherstats=motherstats(~badtracesallPOI,:);
IFvals=IFvals(~badtracesallPOI);
%Cdk2inc=Cdk2inc(~badtracesallPOI); Cdk2low=Cdk2low(~badtracesallPOI);
%%% plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hist(IFvals,2:0.5:12.5); xlim([2 12]);
%hist(IFvals,0:0.03:2.53); xlim([0 2.5]); set(gcf,'color','w','PaperPosition',[0 0 4 4]); saveas(gcf,'h:\Downloads\Fig.jpg');
dscatter(IFvals,tracesgated(:,end-5));
%IFlow=IFvals<1; IFhigh=IFvals>=1;

%hist(POIcalc/framesperhr,10);
%dscatter(tracesgated(:,drugspike-1),tracesgated(:,drugspike+10));
%axis([0.4 1.6 0.4 2.1]);
%set(gcf,'color','w','PaperPosition',[0 0 4 4]);
%saveas(gcf,'h:\Downloads\Fig.jpg');
%%% align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alignPOI=tracestats(:,1); %align by mitosis
%alignPOI=POIcalc;
%alignPOI=ones(size(tracesgated,1),1); %align by first frame
%alignPOI=ones(size(tracesgated,1),1)*drugspike; %align by drugspike
[alignedtime,aligneddata]=aligntraces_4(tracesgated,alignPOI,tracestats,motherstats,daughteroption);
%alignedtime=alignedtime/framesperhr;
alignedtime=(alignedtime)/framesperhr;
%%% plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesmoothoption=1; %1:smoothtraces 0:don't smooth
plottraces_2(aligneddata,alignedtime,xstring,ystring,[ymin ymax],tracesmoothoption,drugspike/framesperhr); %plot traces without categorization
%plottraces_2(aligneddata(IFhigh,:),alignedtime,xstring,ystring,[ymin ymax],tracesmoothoption,drugspike/framesperhr); %plot traces without categorization
%plottraces_2(aligneddata(IFlow,:),alignedtime,xstring,ystring,[ymin ymax],tracesmoothoption,drugspike/framesperhr); %plot traces without categorization
%comparetrends_1(aligneddata,alignedtime,Cdk2inc,Cdk2low,xstring,ystring,[ymin ymax],tracesmoothoption);
%comparetrends_1(aligneddata,alignedtime,IFhigh,IFlow,xstring,ystring,[ymin ymax],tracesmoothoption);
%comparetrends_1(aligneddata,time,trace2high,trace2low,xstring,ystring,ylimits,tracesmoothoption);
%plotpanel(aligneddata,alignedtime,[ymin ymax],tracesmoothoption);

%%% plot IF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if IFoption
%     figure; hold on;
%     IFtime=tracestats(:,3)/framesperhr;
%     %scatter(IFtime(Cdk2inc),IFvals(Cdk2inc),'b','o'); scatter(IFtime(Cdk2low),IFvals(Cdk2low),'r','o','markerfacecolor','r');
%     
%     scatter(IFtime(~isnan(POIcalc)),IFvals(~isnan(POIcalc)),'b','o');
%     scatter(IFtime(isnan(POIcalc)),IFvals(isnan(POIcalc)),'r','o');
%     
%     earlyincorp=isnan(POIcalc) & IFvals>10;
%     plotpanel(aligneddata(earlyincorp,:),alignedtime,[ymin 2],tracesmoothoption);
%     
%     xlabel(xstring);
%     %title('pRb(807/811)'); ylabel('log2(meanRFU)');
%     title('EdU Incorporation'); ylabel('log2(meanRFU)');
%     set(gcf,'color','w','PaperPosition',[0 0 4 7]); %saveas(gcf,'h:\Downloads\FigIF.jpg');
% end
%%% sort POIcalc for heatmap figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sortedPOIcalc,sortidx]=sort(POIcalc);
%[sortedPOIcalc,sortidx]=sort(aligneddata(:,13),'descend'); %sort by value at specific frame
sortedaligneddata=aligneddata(sortidx,:);
%%% generate heatmap (POI must be relative to drugspike/mitosis/etc.) %%%%%
heatmapmin=ymin; heatmapmax=ymax;
%datasorted=saturatepertrace(datasorted);
markoption=1; %1:overlay Cdk2start 2:don't overlay marks
%makeheatmaps_2(sortedaligneddata,sortedPOIcalc,alignedtime,heatmapmax,heatmapmin,markoption,xstring,framesperhr,drugspike);
%%% record number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(aligneddata,1);
fprintf('%0.0f traces\n',numtraces);
end