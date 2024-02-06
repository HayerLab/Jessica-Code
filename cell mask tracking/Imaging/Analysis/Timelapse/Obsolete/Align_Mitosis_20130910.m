function Align_Mitosis(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                traces=main(datadir,shot);
                alltraces=[alltraces;traces];
            end
        end
    end
end
%%% DAPI vs EdU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
framesperhr=5;
numtraces=size(alltraces,1);
maxframe=zeros(numtraces,1);
for i=1:numtraces
    maxframe(i)=find(~isnan(alltraces(i,:)),1,'last');
end
hold on;
for i=1:numtraces
    traceframes=1:maxframe(i);
    plot(traceframes/framesperhr,alltraces(i,traceframes));
end
lastframe=max(maxframe);
xlim([0 lastframe/framesperhr]);
ylim([0 1000]);
xlabel('Time relative to drug spike (hrs)'); ylabel('Sensor (median RFU)');
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
end


function traces=main(datadir,shot)
sensorchannel=6;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are mothers and daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalframes=size(tracedata,2);
samplecells=find(~isnan(tracestats(:,2)));
ismother=zeros(numel(samplecells),1);
for i=1:numel(samplecells)
    s=samplecells(i);
    ismother(i)=ismember(s,tracestats(:,2));
end
samplecells=samplecells(ismother>0);
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numsamples=numel(samplecells);
noise=max(diff(tracedata(samplecells,:,sensorchannel),2,2),[],2)>150;
firstandlast=zeros(numsamples,1);
for i=1:numsamples
    s=samplecells(i);
    firstandlast(i)=tracedata(s,tracestats(s,1),sensorchannel)>100 ...
        && tracedata(s,tracestats(s,3),sensorchannel)>100 ...
        && min(tracedata(s,:,sensorchannel))<50;
end
samplecells=samplecells(~noise & firstandlast);
traces=ones(numel(samplecells),totalframes)*NaN;
for i=1:numel(samplecells)
    s=samplecells(i);
    numframes=tracestats(s,3)-tracestats(s,1)+1;
    traces(i,1:numframes)=tracedata(s,tracestats(s,1):tracestats(s,3),sensorchannel);
end
end
