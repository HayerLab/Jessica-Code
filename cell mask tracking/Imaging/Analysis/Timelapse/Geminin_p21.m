function Geminin_p21(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces_geminin=[];
    alltraces_p21=[];
    alltraces_motherp21=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [traces_geminin,traces_p21,traces_motherp21]=main(datadir,shot);
                alltraces_geminin=[alltraces_geminin;traces_geminin];
                alltraces_p21=[alltraces_p21;traces_p21];
                alltraces_motherp21=[alltraces_motherp21;traces_motherp21];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampletraces_geminin=alltraces_geminin;
sampletraces_p21=alltraces_p21;
sampletraces_motherp21=alltraces_motherp21;
framesperhr=5;
%%% extract geminin risetime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relative=1;
[G1length,badtraces]=getgemininrisetime(sampletraces_geminin,relative);
G1length(badtraces)=[];
sampletraces_p21(badtraces,:)=[];
sampletraces_motherp21(badtraces,:)=[];
numtraces=size(sampletraces_p21,1);
%%% plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keyboard;
hist(G1length);
[~,idx]=sort(G1length);
bluecode=[0 0 1]; blackcode=[0 0 0]; redcode=[1 0 0];
cmap=makecmap(bluecode,blackcode,redcode);
%cidx=round(128*(idx/max(idx))); cidx(cidx<1)=1;
sortedtraces_p21=sampletraces_p21(idx,:);
sortedtraces_motherp21=sampletraces_motherp21(idx,:);
ridx=1:numel(G1length);
cidx=round(128*(ridx/max(ridx))); cidx(cidx<1)=1;
for i=1:numtraces
%for i=1:265
    signal=sortedtraces_p21(i,:);
    if isempty(~isnan(signal))
        continue;
    end
    firstframe=find(~isnan(signal),1,'first');
    lastframe=find(~isnan(signal),1,'last');
    signal=signal(firstframe:lastframe);
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    plot((1:(lastframe-firstframe+1))/framesperhr,signal,'color',cmap(cidx(i),:));
    hold on;
    mothersignal=sortedtraces_motherp21(i,:);
    motherlastframe=find(~isnan(mothersignal),1,'last');
    plot((-1:-1:-motherlastframe)/framesperhr,mothersignal(1:motherlastframe),'color',cmap(cidx(i),:));
    plot([0 0],[min([signal mothersignal]) max([signal mothersignal])],'g');
    G1frame=G1length(idx(i));
    plot(G1frame/framesperhr,signal(G1frame),'ro','markerfacecolor', 'r','markersize',8);
    %axis([-5 G1frame/framesperhr 0 100]);
    %xlim([-5 G1frame/framesperhr]);
end
end


function [traces_geminin,traces_p21,motherp21]=main(datadir,shot)
gemininchannel=6;
p21channel=5;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
numframes=size(tracedata,2);
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(tracestats(:,2)));
samplesignals=tracedata(samplecells,:,gemininchannel);
samplestats=tracestats(samplecells,:);
%%% get cells that are mothers also %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ismother=zeros(numel(samplecells),1);
% for i=1:numel(samplecells)
%     s=samplecells(i);
%     ismother(i)=ismember(s,tracestats(:,2));
% end
% samplecells=samplecells(ismother>0);
% samplesignals=tracedata(samplecells,:,gemininchannel);
% samplestats=tracestats(samplecells,:);
%%% get signals of mother %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
motherp21=ones(numel(samplecells),numframes)*NaN;
mgmax=ones(numel(samplecells),1)*NaN;
for i=1:numel(samplecells)
    mother=samplestats(i,2);
    tempmothp21=tracedata(mother,:,p21channel);
    firstframe=find(~isnan(tempmothp21),1,'first');
    lastframe=find(~isnan(tempmothp21),1,'last');
    motherp21(i,1:(lastframe-firstframe+1))=tempmothp21(lastframe:-1:firstframe);
    mgmax(i)=max(tracedata(mother,:,gemininchannel));
end
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
[newsamplecells,traces_geminin]=gatetraces_geminin(samplecells,samplesignals,samplestats,mgmax,10);
traces_p21=tracedata(newsamplecells,:,p21channel);
newidx=ismember(samplecells,newsamplecells);
motherp21=motherp21(newidx,:);
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