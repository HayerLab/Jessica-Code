function Cdk2_pipdeg(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces_dhb=[];
    alltraces_pip=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [traces_dhb,traces_pip]=main(datadir,shot);
                alltraces_dhb=[alltraces_dhb;traces_dhb];
                alltraces_pip=[alltraces_pip;traces_pip];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampletraces_dhb=alltraces_dhb;
sampletraces_pip=alltraces_pip;
%%% extract Cdk2 features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[minval,badtraces]=getCdk2exitval(sampletraces_dhb);
sampletraces_dhb(badtraces,:)=[];
sampletraces_pip(badtraces,:)=[];
minval(badtraces)=[];
% hist(minval,100); hold on; plot([0.5 0.5],[0 400],'r');
% title('minimum DHB ratio');
% set(gcf,'color','w','PaperPosition',[0 0 4 3]); saveas(gcf,'h:\Downloads\Fig1.jpg');
%%% extract pipdeg features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relative=0;
[degstarts,degends,badtraces]=getdegstartandend_openended(sampletraces_pip,relative);
sampletraces_pip(badtraces,:)=[];
degstarts(badtraces)=[];
degends(badtraces)=[];
sampletraces_dhb(badtraces,:)=[];
minval(badtraces)=[];
%%% Calc G2-phase dhb ratio by min dhb ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keyboard;
G2sample=~isnan(degstarts) & ~isnan(degends);
G2dhbtraces=sampletraces_dhb(G2sample,:);
G2degstarts=degstarts(G2sample);
G2degends=degends(G2sample);
numtraces=numel(G2degstarts);
G2dhbval=ones(numtraces,1)*NaN;
for i=1:numtraces
    G2dhbval(i)=G2dhbtraces(i,G2degends(i));
end
G2minval=minval(G2sample);
%G2inc=G2minval>0.5;
%G2low=G2minval<=0.5;
scatter(G2minval,G2dhbval);

%%% calculate tracelengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=length(minval);
tracelength=ones(numtraces,1)*NaN;
for i=1:numtraces
    goodframes=find(~isnan(sampletraces_pip(i,:)));
    tracelength(i)=range(goodframes);
end
%%% plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
%%%%%% G1 length vs min DHBratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,scatter(degstarts/framesperhr,minval);
axis([0 25 0.35 0.9]);
%%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('G1 duration (hrs)'); ylabel('Minimum DHB ratio (cyto/nuc)');
set(gcf,'color','w','PaperPosition',[0 0 5 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end


function [traces_dhbratio,traces_pip]=main(datadir,shot)
pipchannel=6;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalframes=size(tracedata,2);
samplecells=find(~isnan(tracestats(:,2)));
samplesignals=tracedata(samplecells,:,5);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
samplecells=gatetraces_cdk2(samplecells,samplesignals,samplestats,50);
samplesignals=tracedata(samplecells,:,pipchannel);
samplestats=tracestats(samplecells,:);
[samplecells,traces_pip,samplestats]=gatetraces_pipdegron(samplecells,samplesignals,samplestats);
%%% get data from sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
traces_dhbratio=tracedata(samplecells,:,7)./tracedata(samplecells,:,5);
traces_dhbratio(traces_dhbratio<0)=0;
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