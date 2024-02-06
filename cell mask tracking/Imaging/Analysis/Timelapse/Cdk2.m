function Cdk2(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces_dhb=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                traces_dhb=main(datadir,shot);
                alltraces_dhb=[alltraces_dhb;traces_dhb];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampletraces_dhb=alltraces_dhb;
numtraces=size(sampletraces_dhb,1);
%%% extract Cdk2 features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [minval,badtraces]=getCdk2exitval(sampletraces_dhb);
% minval(badtraces>0)=[]
firstframe=ones(numtraces,1)*NaN;
minval=ones(numtraces,1)*NaN;
for i=1:numtraces
    firstframe(i)=find(~isnan(sampletraces_dhb(i,:)),1);
    minval(i)=sampletraces_dhb(i,firstframe(i)+9);
end
badtraces=minval<0.2 | minval>1.5;
sampletraces_dhb(badtraces,:)=[];
minval(badtraces)=[];
firstframe(badtraces)=[];
%%% plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% histogram of DHB min ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hist(minval,100);
hold on; plot([0.56 0.56],[0 400],'r');
title('DHB ratio 4hrs after mitosis');
set(gcf,'color','w','PaperPosition',[0 0 5 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%%%%%% traces of mitoses at a given time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drugspike=210;
drugsample=find(firstframe>drugspike & firstframe<=(drugspike+4));
drugset_traces=sampletraces_dhb(drugsample,:);
drugset_minval=minval(drugsample);
drugset_firstframe=firstframe(drugsample);
numtraces=length(drugsample);
figure, hold on;
for i=1:numtraces
    if drugset_minval(i)>0.56
        tracecolor='b';
    else
        tracecolor='r';
    end
    lastframe=find(~isnan(drugset_traces(i,:)),1,'last');
    plot(drugset_firstframe(i):lastframe,drugset_traces(i,drugset_firstframe(i):lastframe),tracecolor);
end
plot([drugspike drugspike],[0 2.5],'k'); ylim([0.3 2]);
title('DHB ratio after Cdk4/6i addition');
set(gcf,'color','w','PaperPosition',[0 0 5 3]);
saveas(gcf,'h:\Downloads\Fig.jpg');
end


function traces_dhbratio=main(datadir,shot)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(tracestats(:,2)));
samplesignals=tracedata(samplecells,:,5);
samplestats=tracestats(samplecells,:);
%%% gate: >10frames, positive values, min=zero, normalize, remove noise %%%
samplecells=gatetraces_cdk2(samplecells,samplesignals,samplestats,11);
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