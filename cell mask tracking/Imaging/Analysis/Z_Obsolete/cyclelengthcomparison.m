clear; close all;
load('h:\Downloads\gaptimesdmso.mat','gaptimes');
control = gaptimes;
load('h:\Downloads\gaptimesLY29.mat','gaptimes');
drug = gaptimes;

control = control/5;
drug = drug/5;
control = control(control(:,1)<4,:);
drug = drug(drug(:,1)<4,:);
%plot(control(:,1),control(:,2),'.b','markersize',20);
%hold on;
%plot(drug(:,1),drug(:,2),'.r','markersize',20);
boxplot(axes,[control(:,2);drug(:,2)],[ones(size(control,1),1);2*ones(size(drug,1),1)],'labels',{'DMSO','Harmine'});
set(gcf,'color','w');
%title('Drugspike to Mitosis \n(Spike < 4hrs after previous mitosis)');
saveas(gcf,'h:\Downloads\Fig.jpg');
ttest2(control(:,2),drug(:,2),0.01)
[p,h] = ranksum(control(:,2),drug(:,2),'alpha',0.01)