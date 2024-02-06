clear; close all;
cd ..\Functions; %change directory for function calls
%%%%%%%% define conditions %%%%%%%%%%
conditions = {
    'Media',3,1;
    'DMSO',3,2;
    'Gefitinib',4,1;
    %'U73122',4,3;
    %'Thapsigargin',3,4;
    %'Forskolin',4,5;
    %'Rapamycin',4,6;
    %'U0126',4,7;
    %'LY294002',4,8;
    %'Ruboxistaurin',3,9;
    'Harmine',4,10;
    'Nutlin',4,11;
    %'Ponatinib',4,12;
}
condnum = size(conditions,1);
timeaftermat = [
    1,10;
    %1,20;
    %21,40;
    %41,60;
    %61,80;
    ];
samplenum = size(timeaftermat,1);
bmin = 0.4;
bmax = 1.6;
bstep = (bmax-bmin)/12;
bin = [bmin:bstep:bmax];
figure;
set(gcf,'color','w','PaperPosition',[0 0 8 5]);
for i = 1:condnum
    for j = 1:samplenum
        timeafter = timeaftermat(j,:);
        subplotnum = (i-1)*samplenum + j;
        subaxis(condnum,samplenum,subplotnum,'ML',0.1,'MR',0.1,'MT',0.05,'MB',0.1,'PB',0.05);
        %subplot(condnum,samplenum,subplotnum);
        row = cell2mat(conditions(i,2));
        col = cell2mat(conditions(i,3));
        cd ..\Analysis
        [values,delays] = postdrugmitosis(row,col,timeafter);
        cd ..\Functions
        hist(values,bin);
        low = sum(values<=0.8);
        high = sum(values>0.8);
        lowpercent = low/(low+high);
        if i==1
            start = (timeafter(1)-1)/5;
            finish = timeafter(2)/5;
            title([num2str(start),'-',num2str(finish),'hrs']);
        end
        if j==1
            ylabel(char(conditions(i,1)));
        end
        xlabel([num2str(sprintf('%2.2f',lowpercent*100)),'% below 0.8']);
        xlim([bmin bmax]);
    end
end
saveas(gcf,'h:\Downloads\Fig.jpg');
%close(gcf);
cd ..\Analysis; %return to this directory