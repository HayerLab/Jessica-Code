%function VerifyGating(row,col,site)
%row='C';col='09';site='1';
%row='F';col='02';site='1';
row='B';col='07';site='1';

global rawdir maskdir gatedtraces jitters plotfignum immunoframe nucr
projectpath='H:\Documents\Projects\';
%imagepath='H:\Images\';
imagepath='E:\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130719\';
%experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130905\';
%experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130920\';
experimentpath='2013-06-12_JY_BJ5_CycDFKBP_DHFRp21\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
%shotdir=[num2str(row),'_', num2str(col), '_', num2str(site)];
shotdir=[num2str(row),num2str(col),'\site',num2str(site)];
shotname=[num2str(row),num2str(col),'_site',num2str(site),'_'];
datadir=([projectpath,experimentpath,'Data\']);
rawdir=[imagepath,experimentpath,'Raw\',shotdir,'\',shotname];
maskdir=[imagepath,experimentpath,'Mask\',shotdir,'\',shotname];
immunoframe=0;

%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectmode=1;
selectin=[111]; %C_03_1:101
%selectin=[11,46,69,79,94,143,164,205,230,257,270,310,358];
%selectin=[15,16,50,90,99,169,175,190,198,204,206]; %bin1:bgwin10n
%selectin=[40,59,66,68,82,88,89,98,100];            %bin2
plotsignal2=1;
framesperhr=5;
nucr=4;
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=0; ymax1=1;  %geminin
%ymin1=0.3; ymax1=2.3; %DHB-nuc
%ymin2=0; ymax2=5000;  %geminin
%ymin2=0.3; ymax2=2.3;   %sensor
ymin2=0; ymax2=1;
steps=5;
ystep1=round((ymax1-ymin1)/steps); ystep2=round((ymax2-ymin2)/steps);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8; %default=8 presentation=12
tracewidth=2; %default=2 presentation=4
tracecolor=[.49 1 .63];
drugspike=0/framesperhr; %20130831:154
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% link genealogies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fulltraces=linkgenealogies(tracedata,tracestats);
[totalcells,totalframes,totalsignals]=size(fulltraces);
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstframe=1; lastframe=208; %0719:216 0920:226 CycD-FKBP:208
gatedids=find(~isnan(fulltraces(:,firstframe,1)) & ~isnan(fulltraces(:,lastframe,1)));
gatedtraces=fulltraces(gatedids,:,:);
numgated=size(gatedtraces,1);
%%%%%% signal choice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsig=firstframe:lastframe;
%signal1=gatedtraces(:,xsig,7)./gatedtraces(:,xsig,5);
signal1=gatedtraces(:,xsig,6)./(max(gatedtraces(:,xsig,6),[],2)*ones(1,lastframe-firstframe+1));
signal2=gatedtraces(:,xsig,5)./(max(gatedtraces(:,xsig,5),[],2)*ones(1,lastframe-firstframe+1));
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
xtime=xsig/framesperhr;
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
        %subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
        subaxis(2,3,mod(counter-1,24)+1,'ML',0.04,'MR',0.02,'MT',0.06,'MB',0.06,'SH',0.04,'SV',0.1); %5x4
    end
    set(gcf,'color','w');
    ysig1=signal1(i,:,1);
    if plotsignal2
        ysig2=signal2(i,:,1);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([firstframe/framesperhr lastframe/framesperhr ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1);
        set(hline1,'DisplayName',num2str(i),'color','b','linewidth',tracewidth);
    else
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth);
        axis([firstframe/framesperhr lastframe/framesperhr ymin1 ymax1]);
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike
        line([drugspike drugspike],[ymin1 ymax1],'Color','r','linewidth',tracewidth); 
    end
    %%%%% mark all mitoses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relmitoses=find(gatedtraces(i,xsig,end)==1);
    %absmitoses=xsig(relmitoses);
    plot(xtime(relmitoses),ysig1(relmitoses),'ro','markerfacecolor', 'r','markersize',dotsize);
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([firstframe/framesperhr lastframe/framesperhr ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',tracecolor,'YTick',ymin2:ystep2:ymax2);
        set(hline2,'linewidth',tracewidth);
    end
    title(num2str(i));
    %%% operate image-viewing mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if selectmode
        dcm_obj=datacursormode(gcf);
        datacursormode on;
        response='on';
        while response                                          %get next graph: 'enter'
            response=input('','s');
            figure(gcf);
            if strcmp(response,'1')                             %view angie sensor: 'a','enter'
                delete(findall(gcf,'Type','hggroup'));          %must remove datapoints before changing program
                set(dcm_obj,'Updatefcn',@datatip_image_channel1);
                while ~strcmp(response,'d')                     %exit image mode:'d','enter'
                    response = input('','s');
                end
                delete(findall(gcf,'Type','hggroup'));
                close(gcf+1);                                   %close the image window
            end
            if strcmp(response,'2')                             %view ctd sensor: 'c','enter'
                delete(findall(gcf,'Type','hggroup'));          %must remove datapoints before changing program
                set(dcm_obj,'Updatefcn',@datatip_image_channel2_nomask);
                while ~strcmp(response,'d')                     %exit image mode:'d','enter'
                    response = input('','s');
                end
                delete(findall(gcf,'Type','hggroup'));
                close(gcf+1);                                   %close the image window
            end
            figure(gcf);
        end
    end
end
fprintf([num2str(numgated),'\n']);
%{
title('Sample Trace: Geminin (blue) vs FKBP-CycD1 (green)');
xlabel('Time of Imaging (hrs)'); ylabel('normalized RFU');
set(gcf,'color','w','PaperPosition',[0 0 9 6]); %3x4 or mini
saveas(gcf,'h:\Downloads\Fig.jpg');
close(gcf);
%}