clear; close all;
cd ..\Functions; %change directory for function calls
histview = 0;
trendview = 1;

%%%%%%%% define conditions %%%%%%%%%%
conditions = {
    'Media',3,1;
    'DMSO',3,2;
    'Gefitinib 10uM',4,1;
    %'U73122 0.5uM',3,3;
    %'U73122 5uM',4,3;
    %'Forskolin 1uM',3,5;        %delay
    %'Forskolin 10uM',4,5;        %delay
    %'Thapsigargin 10nM',3,4;
    %'Rap 10nM',3,6;
    %'Rap 100nM',4,6;
    %'U0126 0.1uM',3,7;
    %'U0126 1uM',4,7;
    %'LY29 0.1uM',3,8;
    %'LY29 1uM',4,8;
    %'Ruboxistaurin 1uM',3,9;
    %'Harmine 0.5uM',3,10;
    %'Harmine 5uM',4,10;
    %'Nutlin 1uM',3,11;
    %'Nutlin 10uM',4,11;
    %'Ponatinib 10nM',3,12;
    %'Ponatinib 100nM',4,12;
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels');                %sets screensize units by pixels
screendims = get(0,'ScreenSize');       %get screensize in pixels
screenx = screendims(3);
screeny = screendims(4);

condnum = size(conditions,1);
%timeaftermat = [1,20;21,40;41,60;61,80;];
%timeaftermat = [1,5;6,10;11,15;16,20;21,25;26,30;31,35;36,40;41,45;46,50;51,55;56,60;61,65;66,70;71,75;76,80];
timeaftermat = [1:80]';   %column vector
samplenum = size(timeaftermat,1);
trendval = zeros(condnum,samplenum);

windowaverage = 0;
if windowaverage
    timeaftertitle = cell(1,samplenum);
    for i = 1:samplenum
        start = (timeaftermat(i,1)-1)/5;
        finish = timeaftermat(i,2)/5;
        timeaftertitle{i} = [num2str(start),'-',num2str(finish),'hrs'];
    end
else
    timeaftertitle = timeaftermat/5;
end

if histview
    bmin = 0.4;
    bmax = 1.6;
    bstep = (bmax-bmin)/12;
    bin = [bmin:bstep:bmax];
    histfig = figure;
    %set(histfig,'Position',[0.1 0.5 0.3 0.3]);
    set(gcf,'Position',[round(0.1*screenx) round(0.3*screeny) round(0.4*screenx) round(0.5*screeny)]);
    %set(histfig,'color','w','PaperPosition',[0 0 8 5]); %5x4
    set(histfig,'color','w','PaperPosition',[0 0 8 4]); %3x4 or mini
end

if trendview
    trendfig = figure;
    %set(trendfig,'Position',[0.5 0.5 0.3 0.3]);
    set(gcf,'Position',[round(0.55*screenx) round(0.4*screeny) round(0.3*screenx) round(0.4*screeny)]);
    set(trendfig,'color','w','PaperPosition',[0 0 6 4]);
    if windowaverage
        timeafterordinal = ordinal([1:samplenum],timeaftertitle);
    end
    palette = 'krgcm';
end

for i = 1:condnum
    for j = 1:samplenum
        row = cell2mat(conditions(i,2));
        col = cell2mat(conditions(i,3));
        timeafter = timeaftermat(j,:);
        cd ..\Analysis
        [values,delays] = angiepostdrugmitosis(row,col,timeafter);
        cd ..\Functions
        low = sum(values<=0.8);
        high = sum(values>0.8);
        lowpercent = (low/(low+high))*100;
        highpercent = 100-lowpercent;
        trendval(i,j) = highpercent;
        
        if histview
            figure(histfig);
            subplotnum = (i-1)*samplenum + j;
            %subaxis(condnum,samplenum,subplotnum,'ML',0.1,'MR',0.1,'MT',0.05,'MB',0.1,'PB',0.05); %3x4
            subaxis(condnum,samplenum,subplotnum,'MR',0.05,'MT',0.1,'MB',0.1,'PB',0.05); %mini
            n = hist(values,bin);
            hist(values,bin);
            if i==1
                title(timeaftertitle{j});
            end
            xlabel([num2str(sprintf('%2.2f',lowpercent)),'% below 0.8']);
            xlim([bmin bmax]);
            ylim([0 max(n)]);
            if j==1
                ylabel(char(conditions(i,1)));
            end
        end
    end
    if trendview
        figure(trendfig);
        hold on;
        if windowaverage
            plot(timeafterordinal,trendval(i,:),palette(i),'linewidth',2,'markersize',8,'markerfacecolor',palette(i));
            set(gca,'xtick',1:samplenum);
            set(gca,'xticklabel',timeaftertitle);
        else
            plot(timeaftertitle,trendval(i,:),palette(i),'linewidth',2);
        end
        hold off;
    end
end
if trendview
    figure(trendfig);
    legend(char(conditions(:,1)),'location','east');
    ylabel('cells above 0.8 hDHB ratio (%)');
    xlabel('time after mitosis (hrs)');
end
saveas(gcf,'h:\Downloads\Fig.jpg');
prepforexcel = 0;
if prepforexcel
    conditions(:,1)
    timeaftertitle
    trendval
end
%close(gcf);
cd ..\Analysis; %return to this directory