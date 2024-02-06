function Visualize_CellCycleGating(conditions,datadir)
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    Hoechstvals=[];
    EdUvals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [Hoechstval,EdUval]=main(datadir,shot);
                Hoechstvals=[Hoechstvals;Hoechstval];
                EdUvals=[EdUvals;EdUval];
            end
        end
    end
end
%%% Hoechst histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minx=2000000;maxx=9000000;stepx=(maxx-minx)/100;
% hist(Hoechstvals,minx-stepx:stepx:maxx+stepx);
% xlim([minx maxx]);
% xlabel('Hoechst (sumRFU)'); ylabel('# cells');
% set(gcf,'color','w','PaperPosition',[0 0 4 4]);
% saveas(gcf,'h:\Downloads\Fig.jpg');
%%% EdU histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minx=13;maxx=24;stepx=(maxx-minx)/100;
% hist(EdUvals,minx-stepx:stepx:maxx+stepx);
% xlim([minx maxx]);
% xlabel('EdU log2(sumRFU)'); ylabel('# cells');
% set(gcf,'color','w','PaperPosition',[0 0 2 2]);
% saveas(gcf,'h:\Downloads\Fig.jpg');
%%% DAPI vs EdU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;dscatter(Hoechstvals,EdUvals);
%hist(Hoechstvals,100); xlim([0 600000]);
xlabel('Hoechst (sumRFU)');
ylabel('EdU log2(sumRFU)');
%axis([25000 350000 11 22]);
axis([2000000 9000000 13 24]);

%xlabel('Hoechst (sumRFU)');
%ylabel('#cells');
%axis([12 22 5 10]);

hold on;
%%% DMSO
%rectangle('Position',[2500000,14,2000000,2],'EdgeColor','r','linewidth',2);    %G1
%rectangle('Position',[3000000,21.5,4000000,1.5],'EdgeColor','r','linewidth',2);   %S
rectangle('Position',[6000000,15,2000000,1.5],'EdgeColor','r','linewidth',2);    %G2

%%% 4/6i
% rectangle('Position',[2750000,14,2000000,2],'EdgeColor','r','linewidth',2);    %G1
% rectangle('Position',[3000000,21.5,4000000,1.5],'EdgeColor','r','linewidth',2);   %S
% rectangle('Position',[6250000,15,2000000,1.5],'EdgeColor','r','linewidth',2);    %G2

%%% 2i or (2i+4/6i)
% rectangle('Position',[2750000,13.5,2000000,2],'EdgeColor','r','linewidth',2);    %G1
% rectangle('Position',[3500000,20,3500000,2],'EdgeColor','r','linewidth',2);   %S
% rectangle('Position',[6250000,14.5,2500000,1.5],'EdgeColor','r','linewidth',2);    %G2

%%% All
% rectangle('Position',[2500000,13.5,2250000,2],'EdgeColor','r','linewidth',2);    %G1
% rectangle('Position',[3000000,20,4000000,3],'EdgeColor','r','linewidth',2);   %S
% rectangle('Position',[6000000,14.5,2500000,1.5],'EdgeColor','r','linewidth',2);    %G2


set(gcf,'color','w','PaperPosition',[0 0 4 4]); %default 6 4
saveas(gcf,'h:\Downloads\Fig.jpg');
%}
end


function [Hoechstval,EdUval]=main(datadir,shot)
%load([datadir,shot,'.mat'],'IFdata');
load([datadir,'IF_',shot,'.mat'],'IFdata');
Hoechstval=IFdata(:,3).*IFdata(:,10);
EdUval=log2(IFdata(:,3).*IFdata(:,12)); %sum
%EdUval=log2(IFdata(:,5));
end
