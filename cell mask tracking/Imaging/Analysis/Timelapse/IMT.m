function IMT(conditions,datadir)
condnum=size(conditions,1);
framesperhr=5;
boxplotdata=[];
offset=[];
bstep=1/framesperhr;
bmin=0;
bmax=20;
bin=bmin:bstep:bmax;
namecount=cell(condnum,1);
colorcode='krgc';
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
                [traces,stats]=main(datadir,shot,i);
                alltraces=[alltraces;traces];
                allstats=[allstats;stats];
            end
        end
    end
    %%% remove cells that mitose too late to be considered %%%%%%%%%%%%%%%%%%%%
    totalframes=size(alltraces,2);
    latestframe=totalframes-bmax*framesperhr+1;
    allstats=allstats(allstats(:,1)<=latestframe,:);
    %%% phase lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    values=(allstats(:,3)-allstats(:,1))/framesperhr;
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
title('IMT comparison');
xlabel('duration (hrs)'); ylabel('cdf (%)');
set(gcf,'color','w','PaperPosition',[0 0 4 5]);
saveas(gcf,'h:\Downloads\Fig_cdfs.jpg');
figure,boxplot(axes,boxplotdata,offset,'labels',namecount);
title('IMT comparison');ylabel('duration (hrs)');
set(gcf,'color','w','PaperPosition',[0 0 4 5]);
saveas(gcf,'h:\Downloads\Fig_boxplots.jpg');
end


function [samplesignals,samplestats]=main(datadir,shot,control)
sensorchannel=6;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(tracestats(:,2)));
%%% get cells that are mothers also %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ismother=zeros(numel(samplecells),1);
for i=1:numel(samplecells)
    s=samplecells(i);
    ismother(i)=ismember(s,tracestats(:,2));
end
samplecells=samplecells(ismother>0);
%%% get data from sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplesignals=tracedata(samplecells,:,sensorchannel);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
if control~=1
    [~,samplesignals,samplestats]=gatetraces(samplecells,samplesignals,samplestats);
end
end
