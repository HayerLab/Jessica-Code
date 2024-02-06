function Cdk2(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces_dhbratio=[];
    alltraces_dist=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [traces_dhbratio,traces_dist]=main(datadir,shot);
                alltraces_dhbratio=[alltraces_dhbratio;traces_dhbratio];
                alltraces_dist=[alltraces_dist;traces_dist];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampletraces_dhbratio=alltraces_dhbratio;
sampletraces_dist=alltraces_dist;
%%% extract Cdk2 features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[minval,badtraces]=getCdk2exitval_dist(sampletraces_dhbratio,sampletraces_dist);
sampletraces_dhbratio(badtraces>0,:)=[];
sampletraces_dist(badtraces>0,:)=[];
minval(badtraces>0)=[];
keyboard;
hist(minval,100);
xlabel('minimum DHB ratio');
set(gcf,'color','w','PaperPosition',[0 0 5 5]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%%% extract pipdeg features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[degstarts,degends,badtraces]=getdegstartandend(sampletraces_dist);
sampletraces_dist(badtraces>0,:)=[];
degstarts(badtraces>0)=[];
degends(badtraces>0)=[];
sampletraces_dhbratio(badtraces>0)=[];
minval(badtraces>0)=[];
risetime(badtraces>0)=[];
riseslope(badtraces>0)=[];
%%% calculate additional features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(sampletraces_dist,1);
traceends=ones(numtraces,1)*NaN;
for i=1:numtraces
    traceends(i)=find(~isnan(sampletraces_dist(i,:)),1,'last');
end
slength=degends-degstarts;
g2length=traceends-degends;
%%% plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keyboard;
framesperhr=5;
scatter(degstarts/framesperhr,minval);
axis([0 20 0.35 1]);
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
xlabel('G1 duration (hrs)'); ylabel('Minimum DHB ratio (cyto/nuc)');
set(gcf,'color','w','PaperPosition',[0 0 5 5]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end


function [traces_dhbratio,traces_dist]=main(datadir,shot)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(tracestats(:,2)));
samplesignals=tracedata(samplecells,:,5);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
samplecells=gatetraces_cdk2(samplecells,samplesignals,samplestats,50);
%%% get data from sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces_dhbratio=tracedata(samplecells,:,7)./tracedata(samplecells,:,5);
traces_dhbratio(traces_dhbratio<0)=0;
x=tracedata(samplecells,:,1); y=tracedata(samplecells,:,2);
[numtraces,numframes]=size(x);
traces_dist=ones(numtraces,numframes)*NaN;
nucr=12;
winrad=2.5*nucr;
for i=1:numframes
    otherx=tracedata(:,i,1); otherx=otherx(~isnan(otherx));
    othery=tracedata(:,i,2); othery=othery(~isnan(othery));
    for j=1:numtraces
        xdiff=abs(otherx-x(j,i)); ydiff=abs(othery-y(j,i));
        traces_dist(j,i)=isempty(find(xdiff<winrad & xdiff >0 & ydiff<winrad & ydiff>0,1));
    end
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