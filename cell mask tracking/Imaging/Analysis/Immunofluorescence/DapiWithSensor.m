function DapiWithSensor(conditions,datadir,sensor)
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    samplesize=numel(rowmat)*numel(colmat)*numel(sitemat);
    numcells=zeros(samplesize,1);
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                main(datadir,shot);
            end
        end
    end
end
end

function main(shot)
datadir='H:\Documents\Projects\2013-06-07_p21_cy2_deletions\Experiment_20130719\Data\';
shot='B_07_1';
sensorchannel=6;
IFchannel=9;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
load([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
[totalcells,totalframes,totalsignals]=size(tracedata);
tracedataIF=[tracedata,ones(totalcells,1,totalsignals)*NaN];
totalframes=totalframes+1;
tracedataIF(:,totalframes,:)=IFdata(:,1:totalsignals);
totalframes=size(tracedataIF,2);
staincells=find(~isnan(tracedataIF(:,totalframes,1)));
numstain=numel(staincells);
tracestats(staincells,3)=totalframes;
%%% find traces where sensor is present at any point %%%%%%%%%%%%%%%%%%%%%%
samplecells=find(tracestats(:,3)==totalframes & max(tracedataIF(:,:,sensorchannel),[],2)>=50);
IFval=IFdata(samplecells,IFchannel);
sensorval=IFdata(samplecells,sensorchannel);
dscatter(sensorval,IFval);
%%% find traces where sensor drops before movie ends %%%%%%%%%%%%%%%%%%%%%%
filterlength=10;

for i=1:numel(samplecells)
    if tracedataIF(samplecells(i),tracestats(samplecells(i)),sensorchannel)<50
        samplecells(i)=NaN;
    end
end
samplecells=samplecells(~isnan(samplecells));
numsample=numel(samplecells);

lowfilter = zeros(numsample,totalframes);
droppoi = zeros(numsample,1);
for i=1:numsample
    cellid=staincells(i);
    signal = tracedataIF(cellid,:,sensorchannel);
    reversesignal = signal(totalframes:-1:1);
    rs_slope = getslope_forward(reversesignal,1:filterlength);
    rs_depth = -reversesignal;
    rs_time = 1:totalframes;
    reversefilter = smooth(rs_slope*0.25+rs_depth-rs_time);
    lowfilter(i,:) = reversefilter(totalframes:-1:1);
    droppoi(i) = find(lowfilter(i,:)==max(lowfilter(i,:)),1,'last');
end
timesincedrop = totalframes-droppoi;
%{
scatter(timesincedrop(lowset)/5,EdU(lowset),40);
xlabel('Time since signal drop (hr)');
ylabel('Cyclin D1 Immunofluorescence');
axis([0 15 0 3.5]);
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
%{
xmin=min(timesincedrop(lowset)); xmax=max(timesincedrop(lowset));
x=xmin:xmax;
y=zeros(1,length(x));
e=zeros(1,length(x));
for i=1:length(x)
    index=find(timesincedrop(lowset)==x(i));
    y(i)=mean(EdU(lowset(index)));
    e(i)=std(EdU(lowset(index)));
end
x=x/5;
errorbar(x,y,e,'o','markeredgecolor','k','markerfacecolor',[.49 1 .63],'linewidth',1.5);
xlabel('Time since signal drop (hr)');
ylabel('EdU Incorporation');
axis([0 15 0 3.5]);
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 1;
showsignal2 = 1;
if panelvisualize
    sample=1:samplesize;
    for cc=1:length(sample)
        i=sample(cc);
        figure(ceil(cc/20));             %only 24 plots per figure
        set(gcf,'color','w');
        set(gcf,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]);
        subaxis(4,5,mod(cc-1,20)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
        sensortrace = sensor(i,:);
        if showsignal2==1
            signal2trace = lowfilter(i,:);
            [haxes,hline1,hline2] = plotyy(xtime,sensortrace,xtime,signal2trace);
            axes(haxes(1));
        else
            hline1=plot(xtime,sensortrace);
        end
        hold on;
        axis([0 numframes y1min y1max]);
        set(gca,'YTick',y1min:y1step:y1max);
        set(hline1,'color','b','linewidth',signalwidth);
        title(['Cell ',num2str(i)]);
        
        %scatter(numframes,EdU(i),120,'ro','markerfacecolor','r');
        %scatter(lowpoi(i),sensortrace(lowpoi(i)),80,'ro','markerfacecolor','r');
        scatter(highpoi(i),sensortrace(highpoi(i)),80,'ro','markerfacecolor','r');
        
        if showsignal2==1
            axes(haxes(2)); hold on;
            axis([0 numframes y2min y2max]);
            set(gca,'YTick',y2min:y2step:y2max);
            set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
            %line([xtime(1) xtime(end)],[0 0],'Color','k');
            %}
        end
    end
end
end
%{
scatter(EdU,sensorcalc);
xlabel('EdU Incorporation');
ylabel('p21-sensor');
axis([-0.5 3.5 -0.1 2]);
set(gcf,'color','w','PaperPosition',[0 0 6 9]); %3x4 or mini
saveas(gcf,'h:\Downloads\Fig.jpg');
close(gcf);
%}