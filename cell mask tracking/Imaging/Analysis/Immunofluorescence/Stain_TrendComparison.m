function Stain_TrendComparison(conditions,datadir)
coloroption=1;
colorcharmat={'r','b'};
colorcodemat={[1 0 0],[0 0 1]};
colorchar=colorcharmat{coloroption};
colorcode=colorcodemat{coloroption};

condnum=size(conditions,1);
sigmean=ones(1,condnum)*NaN;
sigstd=ones(1,condnum)*NaN;
sigtime=0:2:16; sigtime2=[sigtime fliplr(sigtime)];
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    Abvals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                Abval=main(datadir,shot);
                Abvals=[Abvals;Abval];
            end
        end
    end
    sigmean(i)=nanmean(Abvals);
    sigstd(i)=nanstd(Abvals);
end
hold on;
color1=colorcode + .7*(1-colorcode);
plot(sigtime,sigmean,colorchar,'linewidth',1.5,'Marker','o','Markersize',4,'Markerfacecolor',colorchar);
fill(sigtime2,[sigmean-sigstd fliplr(sigmean+sigstd)],color1,'edgecolor',color1,'FaceAlpha', 0.4);
titlestring='Cyclin D1';
descstring='IF log2(meanRFU)';
title(titlestring);xlabel('time since serum-release (hr)'); ylabel(descstring);
xlim([min(sigtime) max(sigtime)]); ylim([5 6.2]); %20131209: Goat-pRb(780): 5.4 7
set(gca,'xtick',sigtime);
set(gcf,'color','w','PaperPosition',[0 0 3 3]); %6 4.5
saveas(gcf,'h:\Downloads\Fig.jpg');
end


function Abval=main(datadir,shot)
%1:X 2:Y 3:Area 4:Int(Hoechst) 5:med(pRb) 6:med(p21) 7:med(CycD) 8-10:mean 11-13:max
load([datadir,shot,'.mat'],'IFdata');
Hoechstval=IFdata(:,4);
Abval=log2(IFdata(:,10));
%%% Gate by Hoechst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1cells=Hoechstval<275000; %cycDp21SerumRelease: 300000
Abval=Abval(G1cells);
end
