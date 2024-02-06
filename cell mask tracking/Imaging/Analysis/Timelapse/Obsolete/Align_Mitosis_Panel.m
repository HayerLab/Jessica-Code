function Align_Mitosis_Panel(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces=[];
    allextras=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [traces,extras]=main(datadir,shot);
                alltraces=[alltraces;traces];
                allextras=[allextras;extras];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
sampletraces=alltraces;
sampleextras=allextras;
%%% gate by trace length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces,1);
maxframe=zeros(numtraces,1);
for i=1:numtraces
    maxframe(i)=find(~isnan(sampletraces(i,:)),1,'last');
end
lastframe=prctile(maxframe,50);
shorttraces=maxframe<=lastframe;
sampletraces=sampletraces(shorttraces,:);
sampleextras=sampleextras(shorttraces,:);
maxframe=maxframe(shorttraces);
%%% fill in remainder with zeros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces,1);
for i=1:numtraces
    if maxframe(i)>lastframe
        maxframe(i)=lastframe;
    end
    sampletraces(i,maxframe(i)+1:lastframe)=0;
    sampleextras(i,maxframe(i)+1:lastframe)=0;
end
%%% choose signals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% signal choice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsig=1:lastframe;
xtime=xsig/framesperhr;
signal1=sampleextras(:,xsig);
signal2=sampletraces(:,xsig);
plotsignal2=1;
%%% plot all traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces,1);
selection=1:numtraces;
tracewidth=2;
tracecolor=[.49 1 .63];
for counter=1:length(selection)
    i=selection(counter);
    figure(ceil(counter/24));
    subaxis(4,6,mod(counter-1,24)+1,'ML',0.02,'MR',0.01,'MT',0.03,'MB',0.03,'SH',0.02); %5x4
    set(gcf,'color','w');
    ysig1=signal1(i,:);
    if plotsignal2
        ysig2=signal2(i,:);
        [haxes,hline1,hline2]=plotyy(xtime,ysig1,xtime,ysig2);
        axes(haxes(1));
        %axis([firstframe/framesperhr lastframe/framesperhr ymin1 ymax1]);
        %set(gca,'Box','off','YAxisLocation','left','YTick',ymin1:ystep1:ymax1);
        set(gca,'Box','off','YAxisLocation','left');
        set(hline1,'DisplayName',num2str(i),'color','b','linewidth',tracewidth);
    else
        line(xtime,ysig1,'DisplayName',num2str(i),'linewidth',tracewidth);
        %axis([firstframe/framesperhr lastframe/framesperhr ymin1 ymax1]);
    end
    hold on;
    %%% Adjust signal2 settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotsignal2
        axes(haxes(2));
        %axis([firstframe/framesperhr lastframe/framesperhr ymin2 ymax2]);
        %set(gca,'YAxisLocation','right','YColor',tracecolor,'YTick',ymin2:ystep2:ymax2);
        set(gca,'YAxisLocation','right','YColor',tracecolor);
        set(hline2,'linewidth',tracewidth);
    end
    title(num2str(i));
end
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%xlabel('Time relative to Mitosis (hrs)'); ylabel('Individual traces (sorted by degradation time)');
%set(gcf,'color','w','PaperPosition',[0 0 4 7]);
%saveas(gcf,'h:\Downloads\Fig.jpg');
end


function [traces,extras]=main(datadir,shot)
sensorchannel=5;
extrachannel=6;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalframes=size(tracedata,2);
samplecells=find(~isnan(tracestats(:,2)));
%%% get cells that are also mothers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ismother=zeros(numel(samplecells),1);
% for i=1:numel(samplecells)
%     s=samplecells(i);
%     ismother(i)=ismember(s,tracestats(:,2));
% end
% samplecells=samplecells(ismother>0);
%%% get data from sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplesignals=tracedata(samplecells,:,sensorchannel);
extrasignals=tracedata(samplecells,:,extrachannel);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min>=zero, normalize, remove noise %%
[~,samplesignals,samplestats,extrasignals]=gatetraces_extrachannel(samplecells,samplesignals,samplestats,extrasignals);
%%% align traces to mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(samplesignals,1);
traces=ones(numtraces,totalframes)*NaN;
extras=ones(numtraces,totalframes)*NaN;
for i=1:numtraces
    frames=samplestats(i,1):samplestats(i,3);
    traces(i,1:numel(frames))=samplesignals(i,frames);
    extras(i,1:numel(frames))=extrasignals(i,frames);
end
end

%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
sample=523;
traces=sampletraces(sample,:);
starts=degstarts(sample);
ends=degends(sample);
maxf=maxframe(sample);
samplesize=size(traces,1);
for i=1:samplesize
    %figure(ceil(i/24)); set(gcf,'color','w');
    %subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    x=1:maxf(i);
    plot(x/framesperhr,traces(i,x),'linewidth',4);
    xlim([1 21.6]);
    hold on;
    plot(starts(i)/framesperhr,traces(i,starts(i)),'ro','markerfacecolor','r','markersize',8);
    plot(ends(i)/framesperhr,traces(i,ends(i)),'go','markerfacecolor','g','markersize',8);
end
xlabel('Time relative to Mitosis (hrs)'); ylabel('sensor level');
set(gcf,'color','w','PaperPosition',[0 0 3 2]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%