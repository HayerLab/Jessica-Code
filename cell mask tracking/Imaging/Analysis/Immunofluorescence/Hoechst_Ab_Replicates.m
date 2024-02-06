function Hoechst_EdU_Ab(conditions,datadir)
figure;
for rep=0:1
    condnum=size(conditions,1);
    sigmean=ones(1,condnum)*NaN;
    sigstd=ones(1,condnum)*NaN;
    sigtime=0:11; sigtime2=[sigtime fliplr(sigtime)];
    for i=1:condnum
        rowmat=cell2mat(conditions(i,2))+rep*8;
        colmat=cell2mat(conditions(i,3));
        sitemat=cell2mat(conditions(i,4));
        Ab1vals=[];
        Ab2vals=[];
        Ab3vals=[];
        cc=0;
        for row=rowmat
            for col=colmat
                for site=sitemat
                    cc=cc+1;
                    %shot=wellnum2str(row,col,site);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    [Ab1val,Ab2val,Ab3val]=main(datadir,shot);
                    Ab1vals=[Ab1vals;Ab1val];
                    Ab2vals=[Ab2vals;Ab2val];
                    Ab3vals=[Ab3vals;Ab3val];
                end
            end
        end
        values1=log2(Ab2vals);

        values2=Ab2vals;
        values3=Ab3vals;
        %hist(values2,-5:0.1:10);
        %hist(values1,100);
        %dscatter(values2,values3);

        sigmean(i)=nanmean(values1);
        sigstd(i)=nanstd(values1);

    end
    if rep==0
        s1='r';
        c1=[1 0 0];
    else
        s1='b';
        c1=[0 0 1];
    end
    color1=c1 + .7*(1-c1);
    hold on;
    plot(sigtime,sigmean,s1,'linewidth',1.5);
    %xlim([0 11]); ylim([0 5]);
    fill(sigtime2,[sigmean-sigstd fliplr(sigmean+sigstd)],color1,'edgecolor',color1,'FaceAlpha', 0.4);
end
titlestring='Rabbit pRb(807/811)';
descstring='log2(meanRFU)';
title(titlestring);xlabel('time since serum-release (hr)'); ylabel(descstring);
xlim([0 11]);
set(gcf,'color','w','PaperPosition',[0 0 9 6]); %6 4.5
saveas(gcf,'h:\Downloads\Fig.jpg');
end


function [Ab1,Ab2,Ab3]=main(datadir,shot)
%1:X 2:Y 3:Area 4:Int(Hoechst) 5:med(pRb) 6:med(p21) 7:med(CycD) 8-10:mean 11-13:max
load([datadir,shot,'.mat'],'IFdata');
Hoechstval=IFdata(:,4);
Ab1val=IFdata(:,5); %Ab1val=log2(IFdata(:,5));
Ab2val=IFdata(:,6);
Ab3val=IFdata(:,7);
%Ab3val=IFdata(:,7)./IFdata(:,6); Ab3val(isinf(Ab3val))=NaN;
%%% Gate by Hoechst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1cells=find(Hoechstval<250000); %cycDp21SerumRelease: 300000
Ab1=Ab1val(G1cells);
Ab2=Ab2val(G1cells);
Ab3=Ab3val(G1cells);
end
