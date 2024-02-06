%function VerifyGating(row,col,site)
%row='D';col='05';site='4';
row='B';col='03';site='1';

global rawdir maskdir gatedtraces jitters plotfignum immunoframe nucr
projectpath='H:\Documents\Projects\';
imagepath='H:\Images\';
%experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
experimentpath = '2013-08-05_p21_cy1_deletions\Experiment_20130831\';
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
datadir=([projectpath,experimentpath,'Data\']);
rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
maskdir=[imagepath,experimentpath,'Mask\',shot,'\'];
immunoframe=0;
%%%%%% signal choice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal1=5;
signal2=6;
%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectmode=0;
selectin=[];
plotsignal2=1;
framesperhr=5;
nucr=8;
%%%%%% Axis Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ymin1=-0.2; ymax1=1.3;  %geminin
%ymin1=0; ymax1=1000;
ymin2=-0.2; ymax2=1.3;    %mass
steps=5;
ystep1=round((ymax1-ymin1)/steps); ystep2=round((ymax2-ymin2)/steps);
%%%%%% Viewer Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dotsize=8;
tracewidth=2;
tracecolor=[.49 1 .63];
drugspike=0; %20130831:154
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% link genealogies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fulltraces=linkgenealogies(tracedata,tracestats);
[totalcells,totalframes,totalsignals]=size(fulltraces);
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
firstframe=1; lastframe=totalframes; %totalframes
gatedids=find(~isnan(fulltraces(:,firstframe,1)) & ~isnan(fulltraces(:,lastframe,1)));
gatedtraces=fulltraces(gatedids,:,:);
numgated=size(gatedtraces,1);
%%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels'); screendims=get(0,'ScreenSize'); screenx=screendims(3); screeny=screendims(4);
xsig=firstframe:lastframe;
xtime=xsig/framesperhr;
selection=1:numgated;
if ~isempty(selectin)
    selection=selectin;
end
for counter=selection
    if selectmode
        clf;
        set(gcf,'Position',[round(0.6*screenx) round(0.05*screeny) round(0.38*screenx) round(0.38*screeny)]);
        plotfignum=gcf;
        set(gca,'Position',[0.03 0.1 0.95 0.8]);
    else
        figure(ceil(counter/24));
        subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    end
    set(gcf,'color','w');
    ysig1=gatedtraces(counter,xsig,signal1);
    if plotsignal2
        ysig2=gatedtraces(counter,xsig,signal2);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        axis([firstframe/framesperhr lastframe/framesperhr ymin1 ymax1]);
        set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1);
        set(hline1,'DisplayName',num2str(counter),'color','b','linewidth',tracewidth);
    else
        line(xtime,ysig1,'DisplayName',num2str(counter),'linewidth',tracewidth);
        axis([firstframe/framesperhr lastframe/framesperhr ymin1 ymax1]);
    end
    hold on;
    %%%%% mark drug spike and title %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if drugspike
        line([drugspike drugspike],[ymin1 ymax1],'Color','r'); 
    end
    title(num2str(counter));
    %%%%% mark all mitoses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relmitoses=find(gatedtraces(counter,xsig,end)==1);
    absmitoses=xsig(relmitoses);
    plot(xtime(absmitoses),ysig1(relmitoses),'ro','markerfacecolor', 'r','markersize',dotsize);
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        axis([firstframe/framesperhr lastframe/framesperhr ymin2 ymax2]);
        set(gca,'YAxisLocation','right','YColor',tracecolor,'YTick',ymin2:ystep2:ymax2);
        set(hline2,'color',tracecolor,'linewidth',tracewidth);
    end
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
                set(dcm_obj,'Updatefcn',@datatip_image_signal2);
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
%set(gcf,'color','w','PaperPosition',[0 0 15 10]); %3x4 or mini
%saveas(gcf,'h:\Downloads\Fig.jpg');
%close(gcf);