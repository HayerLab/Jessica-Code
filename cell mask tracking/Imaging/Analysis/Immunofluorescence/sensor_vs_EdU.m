function sensor_vs_EdU(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    sensorvals=[];
    EdUvals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [sensorval,EdUval]=main(datadir,shot);
                sensorvals=[sensorvals;sensorval];
                EdUvals=[EdUvals;EdUval];
            end
        end
    end
end
%%% DAPI vs EdU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
dscatter(log(sensorvals),log(EdUvals));
xlabel('sensor (log(total RFU))');
ylabel('EdU Incorporation (log(max RFU))');
%axis([8.5 13.5 4.5 9.5]); %medEdU
axis([8.5 13.5 5.5 10]); %maxEdU
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
end


function [sensorval,EdUval]=main(datadir,shot)
sensorchannel=6;
signalchannel=6;
EdUchannel=10; %10=max
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
load([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
[totalcells,totalframes,totalsignals]=size(tracedata);
tracedataIF=[tracedata,ones(totalcells,1,totalsignals)*NaN];
totalframes=totalframes+1;
tracedataIF(:,totalframes,:)=IFdata(:,1:totalsignals);
staincells=~isnan(tracedataIF(:,totalframes,1));
tracestats(staincells,3)=totalframes;
%%% find traces where sensor is valid %%%%%%%%%%%%%%%%%%%%%%
samplecells=find(tracestats(:,3)==totalframes & max(tracedataIF(:,:,sensorchannel),[],2)>=50);
sensorval=IFdata(samplecells,signalchannel).*IFdata(samplecells,3);
EdUval=IFdata(samplecells,EdUchannel);
%EdUval=IFdata(samplecells,EdUchannel).*IFdata(samplecells,3);
end
