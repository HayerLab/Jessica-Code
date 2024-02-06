clear; close all;
load('h:\Downloads\intensitiesdmso.mat','intensities');
control = intensities(:);
load('h:\Downloads\intensitiespkcb.mat','intensities');
drug = intensities(:);
boxplot(axes,[control;drug],[ones(length(control),1);2*ones(length(drug),1)],'labels',{'DMSO','PKCbeta'});
set(gcf,'color','w');
%title('Drugspike to Mitosis \n(Spike < 4hrs after previous mitosis)');
saveas(gcf,'h:\Downloads\Fig.jpg');
ttest2(control,drug,0.05)
[p,h] = ranksum(control,drug,'alpha',0.05)