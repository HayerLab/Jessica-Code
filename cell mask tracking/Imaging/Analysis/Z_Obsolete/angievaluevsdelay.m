clear; close all;
cd ..\Functions; %change directory for function calls
%%%%%%%% define conditions %%%%%%%%%%
conditions = {
    'Media',3,1;
    'DMSO',3,2;
    'Gefitinib',4,1;
    'U73122',4,3;
    'Thapsigargin',3,4;
    'Forskolin',4,5;
    'Rapamycin',4,6;
    'U0126',4,7;
    'LY294002',4,8;
    'Ruboxistaurin',3,9;
    'Harmine',4,10;
    %'Nutlin',4,11;
    'Ponatinib',4,12;
}
condnum = size(conditions,1);
%timeafter = [21,40];
timeafter = 20;
figure;
set(gcf,'color','w','PaperPosition',[0 0 8 5]);
filtcoef = 5;
filtfiltvec = (1/filtcoef)*ones(1,filtcoef);
for i = 1:condnum
    %subaxis(3,4,i,'ML',0.1,'MR',0.1,'MT',0.05,'MB',0.1,'PB',0.05);
    subaxis(3,4,i,'PB',0.05);
    %subplot(3,4,i);
    row = cell2mat(conditions(i,2));
    col = cell2mat(conditions(i,3));
    cd ..\Analysis
    [values,delays] = angiepostdrugmitosis(row,col,timeafter);
    cd ..\Functions
    plot(delays,values,'.');
    hold on;
    line(delays,filtfilt(filtfiltvec,1,values),'Color','r','Linewidth',2);
    ylim([0.4 1.6]);
    title(char(conditions(i,1)));
end
saveas(gcf,'h:\Downloads\Fig.jpg');
%close(gcf);
cd ..\Analysis; %return to this directory