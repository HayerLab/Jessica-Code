function ViewTraces(row,col,site)
%row='E';col='03';site='4'; %0715: Same:E_03_1_37, E_03_4_14 Lagging: D_04_2, E_03_4_62
%row='6';col='6';site='2'; %6_6_2 trace 4
%row='D';col='03';site='1'; %C_09_2_10 E_09_1_3
%row='2';col='2';site='1'; %0905: B_06_1_1
row=4;col=8;site=1;

global rawdir maskdir gatedtraces jitters plotfignum immunoframe nucr channel
projectpath='H:\Documents\Projects\';
imagepath='H:\Images\';
%imagepath='E:\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130719\';
%experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130905\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130920\';
%experimentpath='2013-06-12_JY_BJ5_CycDFKBP_DHFRp21\';
%experimentpath='2012-12-11_Sabrina_MCF10A_DHFR-Chy-p21\';
%experimentpath='2013-02-22_Sabrina_MCF10Ap21null_DHFR-Chy-p21\';
%experimentpath='Temporary\ErkCDK2\';
%experimentpath='HeeWon\';
%experimentpath='Kyuho\';
%experimentpath='Sabrina\';
%experimentpath='Steve\';
experimentpath='2014-02-08_H2B-DHB-p21dCy1dK\';

shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
%datadir=([projectpath,experimentpath,'Data\']);
datadir=([projectpath,experimentpath,'Data_chronmerge\']);
separatedirectories=1;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
    maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
else
    rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask\'];
end
immunoframe=0;

%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectmode=0;
selectin=[]; %E_04_4:77 91
plotsignal2=0;
framesperhr=5;
nucr=12;
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ymin1=0; ymax1=1;  %geminin
ymin1=0.3; ymax1=2.5; %DHB
%ymin1=0; ymax1=2000; %abs
%ymin2=0; ymax2=1;
ymin2=1; ymax2=2;
%ymin2=0; ymax2=2000;  %abs
steps=5;
ystep1=round((ymax1-ymin1)/steps); ystep2=round((ymax2-ymin2)/steps);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=2; %default=2 presentation=4
tracecolor=[1 0 0]; %tracecolor=[.49 1 .63];
drugspike=0/framesperhr; %0715:210 0719:149 0831:154 0920:146
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
%%% for purposes of visualizing mis-tracking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reportmistracking(rawdir,nucr,tracedata,genealogy,jitters,tracestats);
%%% link genealogies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fulltraces=linkgenealogies_thin(tracedata,genealogy);
[totalcells,totalframes,totalsignals]=size(fulltraces);
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstframe=1; lastframe=229; %0715:230 0719:216 0831:234 0905:220 0920:226 CycD-FKBP:208 SS1211:230 SS0222:149
xsig=firstframe:lastframe;
gatedids=find(~isnan(fulltraces(:,firstframe,1)) & ~isnan(fulltraces(:,lastframe,1)));
gatedtraces=fulltraces(gatedids,:,:);
%%% additional optional gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% only mitoses %%%%%%%
gatedtraces=gatedtraces(max(gatedtraces(:,xsig,end),[],2)==1,:,:);
samplesignals=gatedtraces(:,:,5);
samplecells=1:size(gatedtraces,1);
%hist(samplesignals(:),0:10:2010); xlim([0 2000]);
%hist(max(samplesignals,[],2),100);
%hist(min(samplesignals,[],2),100);
dimthresh=200;
samplecells=remove_dim_noisy(samplesignals,samplecells,dimthresh); %gate DHBnuc
maxthresh=75; %pipdeg:50 steve:75
minthresh=1000; %pipdeg:25 steve:50
%badtraces=remove_range_noisy(samplesignals,maxthresh,minthresh); samplecells(badtraces)=[];
gatedtraces=gatedtraces(samplecells,:,:);
numgated=size(gatedtraces,1);
%%%%%% signal choice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucarea=gatedtraces(:,xsig,3);
nucdensity=gatedtraces(:,xsig,4)./nucarea;
signal1=gatedtraces(:,xsig,7)./gatedtraces(:,xsig,5);
%signal1=gatedtraces(:,xsig,6);
%signal1=gatedtraces(:,xsig,6).*nucarea;
%signal1=nucdensity;
signal2=gatedtraces(:,xsig,6)./gatedtraces(:,xsig,5);
%signal2=gatedtraces(:,xsig,8);
%signal2=nucarea;
%signal2=nucdensity;
%signal2=gatedtraces(:,xsig,6)./nucarea;
%%% normalize signals
%signal1=signal1./(max(signal1,[],2)*ones(1,lastframe-firstframe+1));
%signal2=signal2./(max(signal2,[],2)*ones(1,lastframe-firstframe+1));
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
xtime=xsig/framesperhr;
if selectmode==0
    drugtime=drugspike;
    drugplot=0;
else
    drugtime=0;
    drugplot=drugspike;
end
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
    else
        figure(ceil(counter/24));
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    end
    set(gcf,'color','w');
    ysig1=smooth(signal1(i,:,1));
    if plotsignal2
        ysig2=smooth(signal2(i,:,1));
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([xtime(1) xtime(end) ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1);
        set(hline1,'DisplayName',num2str(i),'color','b','linewidth',tracewidth);
    else
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth);
        axis([xtime(1) xtime(end) ymin1 ymax1]);
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike
        line([drugplot drugplot],[ymin1 ymax1],'Color','r','linewidth',tracewidth);
        %line([drugtime drugtime],[ymin1 ymax1],'Color','r','linewidth',tracewidth);
    end
    %%%%% mark all mitoses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relmitoses=find(gatedtraces(i,xsig,end)==1);
    %absmitoses=xsig(relmitoses);
    plot(xtime(relmitoses),ysig1(relmitoses),'ro','markerfacecolor', 'r','markersize',dotsize);
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([xtime(1) xtime(end) ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',tracecolor,'YTick',ymin2:ystep2:ymax2);
        set(hline2,'DisplayName',num2str(i),'Color',tracecolor,'linewidth',tracewidth);
    end
    title(num2str(i));
    %%% operate image-viewing mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if selectmode
        dcm_obj=datacursormode(gcf);
        datacursormode on;
        response='on';
        while response                                          %get next graph: 'enter'
            response=input('','s');
            switch response                             %view angie sensor: 'a','enter'
                case '1'
                    channel=1;
                case '2'
                    channel=2;
                case '3'
                    channel=3;
                otherwise
                    response=false;
                    continue;
            end
            figure(gcf);
            delete(findall(gcf,'Type','hggroup'));          %must remove datapoints before changing program
            set(dcm_obj,'Updatefcn',@datatip_image_channelx);
            %set(dcm_obj,'Updatefcn',@datatip_image_channelx_nude);
            while ~strcmp(response,'d')                     %exit image mode:'d','enter'
                response = input('','s');
            end
            delete(findall(gcf,'Type','hggroup'));
            close(gcf+1);                                   %close the image window
            figure(gcf);
        end
    end
end
fprintf([num2str(numgated),'\n']);
%{
title('Sample Trace: DHB ratio (blue) vs DHFR-p21 (green)');
xlabel('Time of Imaging (hrs)'); ylabel('normalized RFU');
set(gcf,'color','w','PaperPosition',[0 0 9 6]); %3x4 or mini
saveas(gcf,'h:\Downloads\Fig.jpg');
close(gcf);
%}