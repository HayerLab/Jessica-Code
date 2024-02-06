function PhaseLengths(conditions,datadir)
condnum=size(conditions,1);
framesperhr=5;
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces=[];
    allstats=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [traces,stats]=main(datadir,shot);
                alltraces=[alltraces;traces];
                allstats=[allstats;stats];
            end
        end
    end
    %%% detect degradation start and end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [degstarts,degends,badtraces]=getdegstartandend_phases(alltraces);
    allstats(badtraces>0,:)=[];
    degstarts(badtraces>0)=[];
    degends(badtraces>0)=[];
    %%% phase lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    notnan=find(~isnan(degends));
    allstats=allstats(notnan,:);
    degstarts=degstarts(notnan,:);
    degends=degends(notnan,:);
    G1lengths=(degstarts-allstats(:,1))/framesperhr;
    Slengths=(degends-degstarts)/framesperhr;
    G2lengths=(allstats(:,3)-degends+1)/framesperhr;
    %%% show data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phasedata=[G1lengths;Slengths;G2lengths];
    offset=[ones(length(G1lengths),1);2*ones(length(Slengths),1);3*ones(length(G2lengths),1)];
    figure,boxplot(axes,phasedata,offset,'labels',{'G1 length','S length','G2 length'});
    %%% display settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ylabel('duration (hrs)');
    set(gcf,'color','w','PaperPosition',[0 0 3 6]);
    saveas(gcf,'h:\Downloads\Fig.jpg');
end
end


function [samplesignals,samplestats]=main(datadir,shot)
sensorchannel=6;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are mothers and daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(tracestats(:,2)));
ismother=zeros(numel(samplecells),1);
for i=1:numel(samplecells)
    s=samplecells(i);
    ismother(i)=ismember(s,tracestats(:,2));
end
samplecells=samplecells(ismother>0);
samplesignals=tracedata(samplecells,:,sensorchannel);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
[~,samplesignals,samplestats]=gatetraces(samplecells,samplesignals,samplestats);
end
