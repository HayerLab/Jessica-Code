function ViewTraces
row=1;col=3;site=1; %7_11_3

global rawdir maskdir tracedata jitters plotfignum immunoframe nucr channel
projectpath='H:\Documents\Projects\';
%imagepath='H:\Images\';
imagepath='E:\';
experimentpath='20131213 R-point CC\20140505 3C G1S TranxTrans\';

shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
datadir=([projectpath,experimentpath,'Data\']);
separatedirectories=1;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\'];
end
immunoframe=0;

%%%%%% Data Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=12;
framesperhr=5;
drugspike=39/framesperhr; %[20 54]
frames=1:103;
%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectmode=0;
selectin=[]; %[82 101 120] [125]
plotsignal2=1;
plotsignal3=0;
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=0.2; ymax1=2;
ymin2=0; ymax2=1; %ymax2=1000;
IFmin1=6; IFmax1=12;
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
%%% for purposes of visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
%reportmistracking_2(rawdir,nucr,tracedata,genealogy,jitters,tracestats);

%%% Select lineage gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motheroption=0; %0:no gating 1:mothers 2:no mothers
daughteroption=1; %0:no gating 1:daughters 2:no daughters
IFoption=0; %0:no IFdata 1:IFdata
[tracedata,tracestats,motherstats,IFdata]=gathertracedata_1(datadir,shot,motheroption,daughteroption,IFoption);

minlengthtrace=40;
minlengthmother=0;
quiescentanalysis=0;
%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucchannel=7; cytochannel=13; %7 13
maxoption=1;
%0:normalize each trace by percentile of entire trace
%1:normalize by first frame of trace (regardless of whether mitosis or not)
maxthresh=100; %threshold above which max of each trace must be
noisethresh=0.2; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
[traces1,badtraces1]=gate_Cdk2_6(tracedata,nucchannel,cytochannel,tracestats,minlengthtrace,maxoption,maxthresh,noisethresh,quiescentanalysis);

if IFoption
    IF1=IFdata(:,9); badtracesIF1=IF1<1; IF1=log2(IF1);
    IF2=IFdata(:,10).*IFdata(:,3); badtracesIF2=IF2<1; IF2=log2(IF2); EdUpos=IF2>16;
end
traces2=tracedata(:,:,8); %mass
%traces3=tracedata(:,:,3); %Area

%%% remove gateout %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
badtraces=badtraces1;
%badtraces=badtraces1 | badtracesIF1 | badtracesIF2 | EdUpos;
traces1=traces1(~badtraces,:);
traces2=traces2(~badtraces,:);
%IF1=IF1(~badtraces,:);
%traces2=traces2(~badtraces,:);
%traces3=traces3(~badtraces,:);
tracedata=tracedata(~badtraces,:,:);
tracestats=tracestats(~badtraces,:);

%%% normalize traces by max in lineage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%0:normalize each trace by max of its own trace
%1:normalize by median of max of all traces if min>max*0.5
%2:normalize each trace by value at first frame after mitosis
%traces1=normalizetraces_3(traces1,tracestats,0);
traces2=normalizetraces_3(traces2,tracestats,0);
%traces3=normalizetraces_3(traces3,tracestats,0);

% lv=traces1(:,end);
% scatter(lv,IF1); axis([ymin1 ymax1 IFmin1 IFmax1]);

numgated=size(tracestats,1);
if numgated>192
    numgated=192;
end
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=2; %default=2 presentation=4
steps=5; ystep1=(ymax1-ymin1)/steps; ystep2=(ymax2-ymin2)/steps; IFstep1=round((IFmax1-IFmin1)/steps);
trace2color=[1 0 0];
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
else
    drugtime=0;
    drugplot=drugspike;
end
xtime=frames/framesperhr;
xtime=xtime-drugtime;
selection=1:numgated;
if ~isempty(selectin)
    selection=selectin;
end
for counter=1:length(selection)
    i=selection(counter);
    if selectmode
        clf;
        set(gcf,'Position',[round(0.6*screenx) round(0.05*screeny) round(0.38*screenx) round(0.38*screeny)]);
        plotfignum=gcf;
        %set(gca,'Position',[0.03 0.1 0.95 0.8]);
    elseif ~isempty(selectin)
        subaxis(1,length(selectin),counter,'ML',0.1,'MR',0.1,'MT',0.1,'MB',0.1,'SH',0.1);
    else
        figure(ceil(counter/24));
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    end
    set(gcf,'color','w');
    ysig1=smoothignorenans(traces1(i,:),3);
    if plotsignal2
        ysig2=smoothignorenans(traces2(i,:),3);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1);
        set(hline1,'DisplayName',num2str(i),'color','b','linewidth',tracewidth);
    elseif IFoption
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime(end),IF1(i));
        axes(haxes(1));
        %hold on; line([xtime(1) xtime(end)],[0.5 0.5],'linestyle','--','color','k');
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1);
        set(hline1,'DisplayName',num2str(i),'color','b','linewidth',tracewidth);
    else
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth);
        axis([xtime(1) xtime(end) ymin1 ymax1]);
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike
        line([drugplot drugplot],[ymin1 ymax1],'Color','k','linewidth',tracewidth,'linestyle','--');
        %line([drugtime drugtime],[ymin1 ymax1],'Color','r','linewidth',tracewidth);
    end
    %%%%% mark mitosis (only the current cell's birth) %%%%%%%%%%%%%%%%%%%%
    relmitosis=tracestats(i,1);
    %plot(xtime(relmitosis),ysig1(relmitosis),'ro','markerfacecolor', 'r','markersize',dotsize);
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([xtime(1) xtime(end) ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',ymin2:ystep2:ymax2);
        set(hline2,'DisplayName',num2str(i),'Color',trace2color,'linewidth',tracewidth);
    end
    if IFoption
        axes(haxes(2));
        hold on;
        scatter(xtime(end),IF1(i),50,'MarkerFaceColor','w','MarkerEdgeColor',trace2color,'linewidth',tracewidth);
        axis([xtime(1) xtime(end) IFmin1 IFmax1]);
        set(gca,'YAxisLocation','right','YColor',trace2color,'YTick',IFmin1:IFstep1:IFmax1);
        %set(hline2,'DisplayName',num2str(i),'MarkerFaceColor','w','MarkerEdgeColor',trace2color,'markersize',12,'linewidth',tracewidth);
    end
    if plotsignal3
        ysig3=smoothignorenans(traces3(i,:),3);
        line(xtime,ysig3,'color','g','DisplayName',num2str(i),'linewidth',tracewidth);
    end
    title(num2str(i));
    %%% operate image-viewing mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if selectmode
        selectmodedisplay;
    end
end
fprintf([num2str(numgated),'\n']);
%{
title('Sample Trace: DHB ratio (blue) vs DHFR-p21 (green)');
xlabel('Time of Imaging (hrs)'); ylabel('normalized RFU');
set(gcf,'color','w','PaperPosition',[0 0 9 3]); %3x4 or mini
saveas(gcf,'h:\Downloads\Fig.jpg');
close(gcf);
%}