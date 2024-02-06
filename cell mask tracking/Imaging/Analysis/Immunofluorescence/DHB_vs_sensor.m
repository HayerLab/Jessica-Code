function DHB_vs_sensor(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    sensorvals=[];
    DHBvals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [sensorval,DHBval]=main(datadir,shot);
                sensorvals=[sensorvals;sensorval];
                DHBvals=[DHBvals;DHBval];
            end
        end
    end
end
%%% DAPI histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
minx=0;maxx=2500000;stepx=(maxx-minx)/100;
hist(DAPIvals,minx:stepx:maxx);
xlim([minx 2000000]);
%}
%%% DAPI vs EdU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
infidx=isinf(DHBvals);
sensorvals=sensorvals(~infidx);
DHBvals=DHBvals(~infidx);
%dscatter(DHBvals,log(sensorvals));
scatter(DHBvals,log(sensorvals),20,'bo','MarkerFaceColor','b');
xlabel('DHB ratio (RFU-cyto/RFU-nuc)');
ylabel('sensor (log(total RFU))');
axis([0.3 2 8 14]);
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
end


function [sensorval,DHBval]=main(datadir,shot)
sensorchannel=6;
signalchannel=11; %10=max 11=mean
nucr=12;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
load([datadir,'sensorIF_',shot,'.mat'],'IFdata','IFjitter');
[totalcells,totalframes,totalsignals]=size(tracedata);
tracedataIF=[tracedata,ones(totalcells,1,totalsignals)*NaN];
totalframes=totalframes+1;
tracedataIF(:,totalframes,:)=IFdata(:,1:totalsignals);
totalframes=size(tracedataIF,2);
staincells=~isnan(tracedataIF(:,totalframes,1));
tracestats(staincells,3)=totalframes;
%%% find traces where sensor is present at any point and not crowded %%%%%%
samplecells=find(tracestats(:,3)==totalframes & max(tracedataIF(:,:,sensorchannel),[],2)>=50);
x=tracedataIF(:,totalframes,1); y=tracedataIF(:,totalframes,2);
winrad=nucr*1;
uncrowded=zeros(numel(samplecells),1);
for c=1:numel(samplecells)
    i=samplecells(c);
    neighbors=find(abs(x-x(i))<winrad & abs(y-y(i))<winrad);
    if numel(neighbors)==1
        uncrowded(c)=1;
    end
end
samplecells=samplecells(uncrowded>0);
%%% calculate values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sensorval=IFdata(samplecells,signalchannel).*IFdata(samplecells,3);
DHBval=IFdata(samplecells,7)./IFdata(samplecells,5);
end
