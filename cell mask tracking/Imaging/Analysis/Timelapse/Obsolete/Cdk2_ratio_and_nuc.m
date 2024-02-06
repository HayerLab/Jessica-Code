function Cdk2_ratio_and_nuc(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces_dhbratio=[];
    alltraces_dhbnuc=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [traces_dhbratio,traces_dhbnuc]=main(datadir,shot);
                alltraces_dhbratio=[alltraces_dhbratio;traces_dhbratio];
                alltraces_dhbnuc=[alltraces_dhbnuc;traces_dhbnuc];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampletraces_dhbratio=alltraces_dhbratio;
sampletraces_dhbnuc=alltraces_dhbnuc;
%%% extract Cdk2 features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[minval,risetime,riseslope,badtraces]=getCdk2features(sampletraces_dhbratio);
sampletraces_dhbratio(badtraces>0,:)=[];
minval(badtraces>0)=[];
risetime(badtraces>0)=[];
riseslope(badtraces>0)=[];
sampletraces_dhbnuc(badtraces>0,:)=[];
%%% extract dhbnuc features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nucval,nucslope,badtraces]=getDHBnucfeatures(sampletraces_dhbnuc,sampletraces_dhbratio);
    
%%% calculate additional features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

traceends=ones(numtraces,1)*NaN;
for i=1:numtraces
    traceends(i)=find(~isnan(sampletraces_dhbnuc(i,:)),1,'last');
end

%%% plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keyboard;
bluecode=[0 0 1]; redcode=[1 0 0]; yellowcode=[1 1 0];
cmap=makecmap(bluecode,yellowcode,redcode);

riseslope_norm=riseslope-prctile(riseslope,5);
riseslope_norm=riseslope_norm/prctile(riseslope_norm,95);
riseslope_color=round(riseslope_norm*128);
riseslope_color(riseslope_color>128)=128;
riseslope_color(riseslope_color<1)=1;

scatter(risetime,minval,10,cmap(riseslope_color,:),'fill');
xmin=prctile(risetime,1); xmax=prctile(risetime,99);
ymin=prctile(minval,1); ymax=prctile(minval,99);
axis([xmin xmax ymin ymax]);

%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('Time relative to Mitosis (hrs)'); ylabel('Individual traces (sorted by degradation time)');
set(gcf,'color','w','PaperPosition',[0 0 4 7]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end


function [traces_dhbratio,traces_dhbnuc]=main(datadir,shot)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalframes=size(tracedata,2);
samplecells=find(~isnan(tracestats(:,2)));
%%% get cells that are also mothers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ismother=zeros(numel(samplecells),1);
for i=1:numel(samplecells)
    s=samplecells(i);
    ismother(i)=ismember(s,tracestats(:,2));
end
samplecells=samplecells(ismother>0);
samplesignals=tracedata(samplecells,:,5);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
samplecells=gatetraces_cdk2(samplecells,samplesignals,samplestats);
samplestats=tracestats(samplecells,:);
%%% get data from sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplesignals_dhbratio=tracedata(samplecells,:,7)./tracedata(samplecells,:,5);
samplesignals_dhbratio(samplesignals_dhbratio<0)=0;
samplesignals_dhbnuc=tracedata(samplecells,:,5);
samplesignals_dhbnuc(samplesignals_dhbnuc<0)=0;
%%% align traces to mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(samplesignals_dhbnuc,1);
traces_dhbratio=ones(numtraces,totalframes)*NaN;
traces_dhbnuc=ones(numtraces,totalframes)*NaN;
for i=1:numtraces
    frames=samplestats(i,1):samplestats(i,3);
    traces_dhbratio(i,1:numel(frames))=samplesignals_dhbratio(i,frames);
    traces_dhbnuc(i,1:numel(frames))=samplesignals_dhbnuc(i,frames)/max(samplesignals_dhbnuc(i,:));
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