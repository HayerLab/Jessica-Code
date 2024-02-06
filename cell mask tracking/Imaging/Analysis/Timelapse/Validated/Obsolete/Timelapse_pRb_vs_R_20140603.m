function Timelapse_pRb_vs_R(conditions,datadir)
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=1; %0:no gating 1:daughters 2:no daughters
quiescentanalysis=1;
if quiescentanalysis
    motheroption=2; daughteroption=2;
end
IFoption=1; %0:no IFdata 1:IFdata
minlengthtrace=54; %CDK2fatecall-minimum=15 Geminin-minimum=10 Pipdeg:5
DHBmin=0.3; DHBmax=2.5;
IFmin=0; IFmax=12;
numbins=50;
IFbin=IFmin:(IFmax-IFmin)/numbins:IFmax;

allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);

xstring='Time relative to Mitosis (hrs)';
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
maxthresh=75; %threshold above which max of each trace must be
noisethresh=0.5; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
[traces1,badtraces1]=gate_Cdk2_1(tracedata,nucchannel,cytochannel,tracestats,minlengthtrace,maxthresh,noisethresh,quiescentanalysis);
%[traces1,badtraces1]=gate_Cdk2_2(tracedata,nucchannel,cytochannel,tracestats,minlengthtrace,maxthresh,noisethresh,quiescentanalysis);

%%% transform and gate IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channelIF=9;
IFvals=IFdata(:,channelIF); badtracesIF=IFvals<1; IFvals=log2(IFvals);
%%% combine gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces=badtraces1 | badtracesIF;
%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces1=traces1(~badtraces,:);
IFvals=IFvals(~badtraces);

%lastval=mean(traces1(:,end-2:end),2);
%lastval=mean(traces1(:,46:48),2);
lastval=traces1(:,end);

% figure;
% dscatter(lastval,IFvals);
% axis([DHBmin DHBmax IFmin IFmax]);
% xlabel('CDK2 activity'); ylabel('pRb log2(meanRFU)');
% set(gcf,'color','w','PaperPosition',[0 0 4 4]);
% saveas(gcf,'h:\Downloads\FigScatter.jpg');

% figure;
% allvals=IFvals(lastval<1);
% hist(allvals,numbins);
% xlim([xmin xmax]);
% set(gcf,'color','w','PaperPosition',[0 0 4 2]);
% saveas(gcf,'h:\Downloads\FigAll.fig');

figure;
lowvals=IFvals(lastval>0.65 & lastval<0.7);
hist(lowvals,IFbin);
xlim([IFmin IFmax]);
set(gcf,'color','w','PaperPosition',[0 0 4 2]);
saveas(gcf,'h:\Downloads\FigLow.fig');

figure;
highvals=IFvals(lastval>0.7 & lastval<75); %[0.7 0.8]
hist(highvals,IFbin);
xlim([IFmin IFmax]);
set(gcf,'color','w','PaperPosition',[0 0 4 2]);
saveas(gcf,'h:\Downloads\FigMidLow.fig');

% figure;
% highvals=IFvals(lastval>1.0 & lastval<1.1);
% hist(highvals,IFbin);
% xlim([IFmin IFmax]);
% set(gcf,'color','w','PaperPosition',[0 0 4 2]);
% saveas(gcf,'h:\Downloads\FigHigh.fig');
end