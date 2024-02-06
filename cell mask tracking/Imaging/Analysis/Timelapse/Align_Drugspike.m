function Align_Drugspike(conditions,datadir)
drugspike=149;
prespikewindow=20;
postspikewindow=20;
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltraces=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                traces=main(datadir,shot,drugspike,prespikewindow,postspikewindow);
                alltraces=[alltraces;traces];
            end
        end
    end
end
%%% DAPI vs EdU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
framesperhr=5;
frames=(-prespikewindow:postspikewindow)/framesperhr;
noise=max(diff(alltraces,2,2),[],2)>150;
alltracesclean=alltraces(~noise,:);
hold on;
for i=1:size(alltracesclean,1)
    plot(frames,alltracesclean(i,:));
end
line([0 0],[0 500],'Color','r','linewidth',2);
ylim([0 500]);
xlabel('Time relative to drug spike (hrs)'); ylabel('Sensor (median RFU)');
set(gcf,'color','w','PaperPosition',[0 0 9 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
end


function traces=main(datadir,shot,drugspike,prespikewindow,postspikewindow)
sensorchannel=6;
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%% gate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(tracedata(:,drugspike,sensorchannel)<=50 ...
    & tracedata(:,drugspike-prespikewindow,sensorchannel)>500 ...
    & ~isnan(tracedata(:,drugspike+postspikewindow,1)));
traces=tracedata(samplecells,drugspike-prespikewindow:drugspike+postspikewindow,sensorchannel);
end
