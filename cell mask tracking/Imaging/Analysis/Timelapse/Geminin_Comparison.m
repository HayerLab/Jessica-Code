function Geminin_Comparison(conditions,datadir)
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
    alltraces_geminin=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                traces_geminin=main(datadir,shot,i);
                alltraces_geminin=[alltraces_geminin;traces_geminin];
            end
        end
    end
    %%% remove cells that mitose too late to be considered %%%%%%%%%%%%%%%%%%%%
    sampletraces_geminin=alltraces_geminin;
    %%% extract geminin risetime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    relative=1;
    [G1length,badtraces]=getgemininrisetime(sampletraces_geminin,relative);
    G1length(badtraces)=[];
    %%% phase lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    values=G1length/framesperhr;
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
title('G1 comparison');
xlabel('duration (hrs)'); ylabel('cdf (%)');
set(gcf,'color','w','PaperPosition',[0 0 4 5]);
saveas(gcf,'h:\Downloads\Fig_cdfs.jpg');
figure,boxplot(axes,boxplotdata,offset,'labels',namecount);
title('G1 comparison');ylabel('duration (hrs)');
set(gcf,'color','w','PaperPosition',[0 0 4 5]);
saveas(gcf,'h:\Downloads\Fig_boxplots.jpg');
end


function [samplesignals,samplestats]=main(datadir,shot,control)
gemininchannel=6;
pipchannel=5;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcells=size(tracedata,1);
tracestats=ones(numcells,3)*NaN;
for c=1:size(tracedata,1)
    tracestats(c,1)=find(~isnan(tracedata(c,:,1)),1,'first');
    tracestats(c,3)=find(~isnan(tracedata(c,:,1)),1,'last');
end
tracestats(:,2)=genealogy;
samplecells=find(~isnan(tracestats(:,2)));
samplestats=tracestats(samplecells,:);
%%% gate by proper pip sensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if control~=1
    pipsignals=tracedata(samplecells,:,pipchannel);
    [samplecells,~,samplestats]=gatetraces_pipdegron(samplecells,pipsignals,samplestats);
end
samplesignals=tracedata(samplecells,:,gemininchannel);
%%% get signal max of mother
sigmax=ones(numel(samplecells),1)*NaN;
for i=1:numel(samplecells)
    mother=samplestats(i,2);
    sigmax(i)=max(tracedata(mother,:,gemininchannel));
end
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
[samplecells,samplesignals]=gatetraces_geminin(samplecells,samplesignals,samplestats,sigmax,50);
samplestats=tracestats(samplecells,:);
end
