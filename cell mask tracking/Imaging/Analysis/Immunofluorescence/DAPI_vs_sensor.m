function DAPI_vs_sensor(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    DAPIvals=[];
    sensorvals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [DAPIval,sensorval]=main(datadir,shot);
                DAPIvals=[DAPIvals;DAPIval];
                sensorvals=[sensorvals;sensorval];
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
dscatter(DAPIvals,log(sensorvals));
xlabel('DAPI (total RFU)');
ylabel('sensor (log(total RFU))');
axis([200000 1500000 9 14]);
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
end


function [DAPIval,sensorval]=main(datadir,shot)
sensorchannel=6;
DAPIchannel=8;
signalchannel=11; %10=max 11=mean
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
%%% find traces where sensor is present at any point %%%%%%%%%%%%%%%%%%%%%%
samplecells=find(tracestats(:,3)==totalframes & max(tracedataIF(:,:,sensorchannel),[],2)>=50);
DAPIval=IFdata(samplecells,DAPIchannel).*IFdata(samplecells,3);
%sensorval=IFdata(samplecells,signalchannel);
sensorval=IFdata(samplecells,signalchannel).*IFdata(samplecells,3);
end
