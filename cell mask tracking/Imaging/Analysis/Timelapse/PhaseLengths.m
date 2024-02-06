function PhaseLengths(conditions,datadir)
condnum=size(conditions,1);
framesperhr=5;
boxplotdata=[];
offset=[];
bstep=1/framesperhr;
bmin=0;
bmax=30;
bin=bmin:bstep:bmax;
namecount=cell(condnum,1);
colorcode='krg';
figure;hold on;
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
    [degstarts,degends,badtraces]=getdegstartandend(alltraces);
    allstats(badtraces>0,:)=[];
    degstarts(badtraces>0)=[];
    degends(badtraces>0)=[];
    %%% remove cells that mitose too late to be considered %%%%%%%%%%%%%%%%%%%%
    totalframes=size(alltraces,2);
    latestframe=totalframes-bmax*framesperhr+1;
    samplecells=allstats(:,1)<=latestframe; %G1:allstats(:,1) S:degstarts G2:degends
    allstats=allstats(samplecells,:);
    degstarts=degstarts(samplecells,:);
    degends=degends(samplecells,:);
    %%% phase lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    notnan=find(~isnan(degends));
    allstats=allstats(notnan,:);
    degstarts=degstarts(notnan,:);
    degends=degends(notnan,:);
    G1lengths=(degstarts-allstats(:,1))/framesperhr;
    Slengths=(degends-degstarts)/framesperhr;
    G2lengths=(allstats(:,3)-degends+1)/framesperhr;
    values=G1lengths;
    %%% store boxplot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    boxplotdata=[boxplotdata;values];
    offset=[offset;i*ones(length(values),1)];
    %%% plot cdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_elements = histc(values,bin);
    c_elements = cumsum(n_elements);
    c_elements = 100*c_elements/max(c_elements);
    c_elements(c_elements==max(c_elements))=99;
    stairs(bin,c_elements,colorcode(i),'linewidth',2);
    namecount{i}=[char(conditions(i,1))];
end
xlim([bmin bmax]);
legend(char(namecount(:)),'location','southeast');
title('G1 duration comparison');
xlabel('duration (hrs)'); ylabel('cdf (%)');
set(gcf,'color','w','PaperPosition',[0 0 4 5]);
saveas(gcf,'h:\Downloads\Fig_cdfs.jpg');
figure,boxplot(axes,boxplotdata,offset,'labels',namecount);
title('G1 duration comparison');ylabel('duration (hrs)');
set(gcf,'color','w','PaperPosition',[0 0 4 5]);
saveas(gcf,'h:\Downloads\Fig_boxplots.jpg');
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
[~,samplesignals,samplestats]=gatetraces_pipdegron(samplecells,samplesignals,samplestats);
end
