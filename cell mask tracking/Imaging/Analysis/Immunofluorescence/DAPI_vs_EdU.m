function DAPI_vs_EdU(conditions,datadir,sensor)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    DAPIvals=[];
    EdUvals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                [DAPIval,EdUval]=main(datadir,shot,sensor);
                DAPIvals=[DAPIvals;DAPIval];
                EdUvals=[EdUvals;EdUval];
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
dscatter(DAPIvals,log(EdUvals));
xlabel('DAPI (total RFU)');
ylabel('EdU Incorporation (log(median RFU))');
axis([200000 1400000 4.5 9.5]); %median
axis([200000 1400000 5.5 9.8]); %max
axis([200000 1400000 11 15.5]); %int
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
end


function [DAPIval,EdUval]=main(datadir,shot,sensor)
sensorchannel=6;
DAPIchannel=8;
EdUchannel=9; %9=median 10=max 11=mean
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
load([datadir,'IF_',shot,'.mat'],'IFdata','IFjitter');
[totalcells,totalframes,totalsignals]=size(tracedata);
tracedataIF=[tracedata,ones(totalcells,1,totalsignals)*NaN];
totalframes=totalframes+1;
tracedataIF(:,totalframes,:)=IFdata(:,1:totalsignals);
totalframes=size(tracedataIF,2);
staincells=~isnan(tracedataIF(:,totalframes,1));
tracestats(staincells,3)=totalframes;
%%% find traces where sensor is present at any point %%%%%%%%%%%%%%%%%%%%%%
if sensor
    samplecells=find(tracestats(:,3)==totalframes & max(tracedataIF(:,:,sensorchannel),[],2)>=50);
else
    samplecells=find(tracestats(:,3)==totalframes);
end
DAPIval=IFdata(samplecells,DAPIchannel).*IFdata(samplecells,3);
EdUval=IFdata(samplecells,EdUchannel);
%EdUval=IFdata(samplecells,EdUchannel).*IFdata(samplecells,3);
end
