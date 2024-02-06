function pRbVsCDK2(conditions,datadir)
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=1; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=0;
if quiescentanalysis
    motheroption=2; daughteroption=2;
end
IFoption=1; %0:no IFdata 1:IFdata
drugspike=0; %3C_CDK4iWee1i:14 %CyclingCDK4iWee1iCdc25i:59
minlengthtrace=30; %CDK2fatecall-minimum=15 Geminin-minimum=10 Pipdeg:5
minlengthmother=0; %PipDeg:50
framesperhr=5;
ymin=0.3; ymax=2.5; ystring='CDK2 activity';

allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);

xstring='Time relative to Mitosis (hrs)';
%xstring='Time relative to drugspike (hrs)';
%xstring='Time since serum release (hrs)';
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
nucchannel=6; cytochannel=8; %6 8
maxthresh=1000; %threshold above which max of each trace must be. %50
noisethresh=0.5; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
[traces1,badtraces1]=gate_Cdk2_1(tracedata,nucchannel,cytochannel,tracestats,minlengthtrace,maxthresh,noisethresh,quiescentanalysis);

%%% transform and gate IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFvals=ones(size(badtraces1))*NaN;
if IFoption
    channelIF=9;
    IFvals=IFdata(:,channelIF);
    badtracesIF=IFvals<1; IFvals=log2(IFvals);
    %badtracesIF=IFvals<0;
end
%%% combine gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces=badtraces1 | badtracesIF;
%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces1=traces1(~badtraces,:);
tracestats=tracestats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
IFvals=IFvals(~badtraces);

% dscatter(traces1(:,end),IFvals);
% axis([0.2 1.4 0 12]);
% xlabel('CDK2 activity'); ylabel('pRb log2(meanRFU)');
% set(gcf,'color','w','PaperPosition',[0 0 4 4]);
% saveas(gcf,'h:\Downloads\FigScatter.jpg');

% figure;
% bmin=0; bmax=12;
% bstep=(bmax-bmin)/50;
% bin=bmin:bstep:bmax;
% numbins=numel(bin);
% binfill=[bin fliplr(bin)];
% Abvals=IFvals(traces1(:,end)<0.85);
% pdfvals=histc(Abvals,bin);
% pdfvals=100*pdfvals/sum(pdfvals);
% fill(binfill,[pdfvals;zeros(numbins,1)],'r','edgecolor','r','FaceAlpha', 0.7);
% xlabel('pRb log2(meanRFU)'); ylabel('pdf (%)');
% set(gcf,'color','w','PaperPosition',[0 0 4 4]);
% saveas(gcf,'h:\Downloads\FigHist.jpg');

% figure;
% allvals=IFvals(traces1(:,end)<0.85);
% hist(allvals,25);
% xlim([0 12]);
% set(gcf,'color','w','PaperPosition',[0 0 4 2]);
% saveas(gcf,'h:\Downloads\FigAll.jpg');
% 
% figure;
% lowvals=IFvals(traces1(:,end)<0.5);
% hist(lowvals,25);
% xlim([0 12]);
% set(gcf,'color','w','PaperPosition',[0 0 4 2]);
% saveas(gcf,'h:\Downloads\FigLow.jpg');
% 
% figure;
% highvals=IFvals(traces1(:,end)>0.6 & traces1(:,end)<0.85);
% hist(highvals,25);
% xlim([0 12]);
% set(gcf,'color','w','PaperPosition',[0 0 4 2]);
% saveas(gcf,'h:\Downloads\FigHigh.jpg');

%%% get Points Of Interest (absolute rather than relative to mitosis) %%%%%
[POI_Cdk2,badtraces_Cdk2]=GetPointsOfInterest_2('Cdk2',traces1,tracestats,minlengthtrace); %detect Cdk2start for daughtertraces
POIcalc=POI_Cdk2-tracestats(:,1)+1;
%%% determine POI gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%badtracesallPOI=zeros(size(POIcalc,1),1); %no bad trace gating
badtracesallPOI=badtraces_Cdk2 | POIcalc<40;
%%% gate out badPOI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
POIcalc=POIcalc(~badtracesallPOI);
tracesgated=traces1(~badtracesallPOI,:); tracestats=tracestats(~badtracesallPOI,:); motherstats=motherstats(~badtracesallPOI,:);
IFvals=IFvals(~badtracesallPOI);
%%% plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(POIcalc/framesperhr,10);
%dscatter(tracesgated(:,drugspike-1),tracesgated(:,drugspike+10));
%axis([0.4 1.6 0.4 2.1]);
%set(gcf,'color','w','PaperPosition',[0 0 4 4]);
%saveas(gcf,'h:\Downloads\Fig.jpg');
%%% align traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cdk2inc=Cdk2inc(~badtracesallPOI); Cdk2low=Cdk2low(~badtracesallPOI);
alignPOI=tracestats(:,1); %align by mitosis
[alignedtime,aligneddata]=aligntraces_4(tracesgated,alignPOI,tracestats,motherstats,daughteroption);
alignedtime=alignedtime/framesperhr;
%%% plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesmoothoption=1; %1:smoothtraces 0:don't smooth
%plottraces_2(aligneddata,alignedtime,xstring,ystring,[ymin ymax],tracesmoothoption,drugspike/framesperhr); %plot traces without categorization

%%% plot IF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure; hold on;
% IFtime=tracestats(:,3)/framesperhr;
% scatter(IFtime(~isnan(POIcalc)),IFvals(~isnan(POIcalc)),'b','o');
% scatter(IFtime(isnan(POIcalc)),IFvals(isnan(POIcalc)),'r','o');
% earlyincorp=isnan(POIcalc) & IFvals>10;
% plotpanel(aligneddata(earlyincorp,:),alignedtime,[ymin 2],tracesmoothoption);
% xlabel(xstring);
% title('pRb(807/811)'); ylabel('log2(meanRFU)');
% set(gcf,'color','w','PaperPosition',[0 0 4 7]); %saveas(gcf,'h:\Downloads\FigIF.jpg');
%%% record number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(aligneddata,1);
fprintf('%0.0f traces\n',numtraces);
end