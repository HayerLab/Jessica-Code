clear; close all;
cd ..\Functions; %change directory for function calls
histview = 0;
stairview = 1;

%%%%%%%% define conditions %%%%%%%%%%
%{
conditions = {
    'Media',3,1;
    'DMSO',3,2;
    %'Gefitinib 10uM',4,1;        %full block
    %'U73122 0.5uM',3,3;           %early accelerate
    %'U73122 5uM',4,3;           %early accelerate + late block
    %'Forskolin 1uM',3,5;        %delay
    %'Forskolin 10uM',4,5;        %delay
    %'Thapsigargin 10nM',3,4;     %big delay
    %'Rapamycin 10nM',3,6;        %delay + early delay
    %'Rapamycin 100nM',4,6;        %delay + early delay
    %'U0126 0.1uM',3,7;            %delay + late delay
    %'U0126 1uM',4,7;            %delay + late delay
    %'LY294002 0.1uM',3,8;         %synchronizer
    %'LY294002 1uM',4,8;         %synchronizer
    %'Ruboxistaurin',3,9;    %early delay
    %'Harmine 0.5uM',3,10;         %late block
    %'Harmine 5uM',4,10;         %late block
    'Nutlin 1uM',3,11;          %full block
    %'Nutlin 10uM',4,11;          %full block
    %'Ponatinib 10nM',3,12;       %accelerate
    %'Ponatinib 100nM',4,12;       %accelerate
}
%}
conditions = {
    %'Non-specific siRNA',[3 4 5],1; %6.5-30.5hr
    %'p21 siRNA',[3 4 5],11;         %6.5-30.5hr
    %'Non-specific siRNA',[5 6 7],1; %24-48hr
    %'p21 siRNA',[5 6 7],11;         %24-48hr
    'Non-specific siRNA',[5 6 7],1; %48-72hr
    'p21 siRNA',[5 6 7],10;         %48-72hr    
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels');                %sets screensize units by pixels
screendims = get(0,'ScreenSize');       %get screensize in pixels
screenx = screendims(3);
screeny = screendims(4);

condnum = size(conditions,1);
samplesize = zeros(condnum);
namecount = cell(condnum,1);
bmin = 0;
bmax = 16;
bstep = (bmax-bmin)/(bmax*5);
bin = [bmin:bstep:bmax];

if histview
    histfig = figure;
    set(histfig,'Position',[round(0.1*screenx) round(0.3*screeny) round(0.3*screenx) round(0.5*screeny)]);
    %set(histfig,'color','w','PaperPosition',[0 0 8 5]); %5x4
    %set(histfig,'color','w','PaperPosition',[0 0 8 4]); %3x4
    %set(histfig,'color','w','PaperPosition',[0 0 2 6]); %Cdt 4x1
    set(histfig,'color','w','PaperPosition',[0 0 3 6]); %Cdt 3x1
    %set(histfig,'color','w','PaperPosition',[0 0 3 4]); %Cdt 2x1
end
if stairview
    stairfig = figure;
    set(stairfig,'Position',[round(0.4*screenx) round(0.3*screeny) round(0.3*screenx) round(0.3*screeny)]);
    set(stairfig,'color','w','PaperPosition',[0 0 6 4]); %cumulative stairs
    palette = 'krgc';
end

for i = 1:condnum
    hold on;
    row = cell2mat(conditions(i,2));
    col = cell2mat(conditions(i,3));
    cd ..\Analysis
    [values,noentry,delays] = cdt1postdrugmitosis(row,col);
    samplesize(i) = size(values,1);
    cd ..\Functions
    entry = numel(values);
    values = [values;(bmax-1)*ones(noentry,1)];
    if histview
        figure(histfig);
        %subaxis(3,1,i,'SV',0.05,'MR',0.05,'MT',0.1,'MB',0.1,'PB',0.07); %5x4
        %subaxis(4,1,i,'MR',0.05,'MT',0.1,'MB',0.1,'PB',0.05); %Cdt 4x1
        subaxis(3,1,i,'SV',0.08,'ML',0.2,'MR',0.05); %Cdt 3x1
        %subaxis(2,1,i,'SV',0.08,'ML',0.2,'MR',0.05,'PT',0.1); %Cdt 2x1        
        hist(values,bin);
        noentrypercent = 100*noentry/(noentry+entry);
        xlabel([sprintf('%0.2f',noentrypercent),'% arrested']);
        ylabel('total cells entered into S-phase (%)');
    end
    if stairview
        figure(stairfig);
        n_elements = histc(values,bin);
        c_elements = cumsum(n_elements);
        c_elements = 100*c_elements/max(c_elements);
        c_elements(c_elements==max(c_elements))=99;         %stupid, but can't see top line otherwise
        stairs(bin,c_elements,palette(i),'linewidth',2);
    end
end
if stairview
    figure(stairview);
    for i=1:condnum
        namecount{i} = [char(conditions(i,1)),' (n = ',mat2str(samplesize(i)),')'];
    end
    legend(char(namecount(:)),'location','southeast');
    %legend(char(conditions(:,1)),'location','southeast');
    title('Mitosis-->Cdt1 degradation');
    xlabel('time of Cdt1 degradation after mitosis (hrs)');
    ylabel('total cells entered into S-phase (%)');    
end
saveas(gcf,'h:\Downloads\Fig.jpg');
%close(histfig);
cd ..\Analysis; %return to this directory