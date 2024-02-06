clear; close all;
angiemode = 1;
cdt1mode = 0;
%%%%%%%% define conditions %%%%%%%%%%
conditions = {
    'Media',3,1;
    'DMSO',3,2;
    
    'Gefitinib 10uM',4,1;
    'U73122 0.5uM',3,3;
    'U73122 5uM',4,3;
    'Thapsigargin 10nM',3,4;
    'Forskolin 1uM',3,5;
    'Forskolin 10uM',4,5;
    'Rapamycin 10nM',3,6;
    'Rapamycin 100nM',4,6;
    'U0126 0.1uM',3,7;
    'U0126 1uM',4,7;
    %}
    %{
    'LY294002 0.1uM',3,8;
    'LY294002 1uM',4,8;
    'Ruboxistaurin',3,9;
    'Harmine 0.5uM',3,10;
    'Harmine 5uM',4,10;
    'Nutlin 1uM',3,11;
    %'Nutlin 10uM',4,11;
    'Ponatinib 10nM',3,12;
    'Ponatinib 100nM',4,12;
    %}
}
if angiemode
    delaypoint = 5;
end
condnum = size(conditions,1);
figure;
set(gcf,'color','w','PaperPosition',[0 0 8 5]);
filtcoef = 5;
filtfiltvec = (1/filtcoef)*ones(1,filtcoef);
for i = 1:condnum
    cd ..\Functions;
    %subaxis(3,4,i,'ML',0.1,'MR',0.1,'MT',0.05,'MB',0.1,'PB',0.05);
    subaxis(3,4,i,'PB',0.05);
    %subplot(3,4,i);
    cd ..\Analysis;
    row = cell2mat(conditions(i,2));
    col = cell2mat(conditions(i,3));
    if angiemode
        [values,delays] = angiepostdrugmitosis(row,col);
        values = values(:,delaypoint);
        goodcells = find(values>=0);
        values = values(goodcells);
        delays = delays(goodcells);
        plot(delays,values,'.');
    end
    if cdt1mode
        [values,noentry,delays] = cdt1postdrugmitosis(row,col);
        plot(delays,values,'.');
    end
    hold on;
    line(delays,filtfilt(filtfiltvec,1,values),'Color','r','Linewidth',2);
    %line(delays,tsmovavg(values,'s',10),'Color','r','Linewidth',2);
    maxdelay = max(delays);
    if isempty(maxdelay)
        maxdelay = 20;
    end
    xlim([0 maxdelay]); 
    ylim([0.3 1.5]);
    %xlim([0 6]); ylim([0 20]);
    title(char(conditions(i,1)));
end
saveas(gcf,'h:\Downloads\Fig.jpg');
%close(gcf);