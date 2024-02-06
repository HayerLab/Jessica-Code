function Histogram_Comparison(conditions,datadir)
condnum=size(conditions,1);
%colorcode=distributecolormap(jet,condnum);
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
    Abvals;
    if dataset==20131126
        %maxz=6.5; minz=4.5; midz=5.75; %serumrelease 20131126 Goat780
        maxz=6.5; minz=4.5; midz=5.5; %serumrelease 20131126 Goat807/811
    elseif dataset==20131209
        maxz=7.5; minz=5; midz=6.25; %serumrelease 20131126 (10perc10nucr median)
    end
end
title(char(conditions(i,1)),'FontSize',16);
%xstring='p21 log2(meanRFU)'; ystring='CycD1 log2(meanRFU)';
%xlabel(xstring); ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 6 5]); %6 5 (small panel: 4 3.3)
saveas(gcf,'h:\Downloads\Fig.jpg');
end


function Abval=main(datadir,shot)
load([datadir,shot,'.mat'],'IFdata');
%Hoechstval=IFdata(:,3).*IFdata(:,8);
Hoechstval=IFdata(:,4);
Ab1val=IFdata(:,8); %perc10nucr10:6 org5 med5 mean9 max13
Ab1val(Ab1val<1)=1; Ab1val=log2(Ab1val);
%%% Gate by Hoechst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1cells=Hoechstval<275000; %pRbAbChar 20131217 Bgsub G1:<250000 G2:400000<x<600000
Abval=Ab1val(G1cells);
end