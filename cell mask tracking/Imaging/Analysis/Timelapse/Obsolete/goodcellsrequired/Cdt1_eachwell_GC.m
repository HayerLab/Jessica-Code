clear; close all;
cd ..\Functions; %change directory for function calls
histview = 0;
stairview = 0;
scatterview = 1;

%%%%%%%% define conditions %%%%%%%%%%

conditions = {
    'DMSO All',1:6,10;
    
    %'DMSO 1',1,10;
    %'DMSO 2',2,10;
    %'DMSO 3',3,10;
    %'DMSO 4',4,10;
    %'DMSO 5',5,10;
    %'DMSO 6',6,10;
    %'DMSO 7',7,10;
    %'DMSO 8',8,10;
    %}
    'MK-2206 1uM 1',1,4;
    'MK-2206 1uM 2',2,4;
    'MK-2206 1uM 3',3,4;
    'MK-2206 1uM 4',4,4;
    'MK-2206 1uM 5',5,4;
    %'MK-2206 1uM 6',6,4; %low sample
    'MK-2206 1uM 7',7,4;
    'MK-2206 1uM 8',8,4;
    %}
    %{
    'Harmine 1uM  1',1,8;
    'Harmine 1uM  2',2,8;
    %'Harmine 1uM  3',3,8; %low sample
    'Harmine 1uM  4',4,8;
    'Harmine 1uM  5',5,8;
    %'Harmine 1uM  6',6,8; %low sample
    %'Harmine 1uM  7',7,8; %low sample
    %'Harmine 1uM 8',8,8; %low sample
    %}
    %{
    'ATMi 1x 1',1,9;
    'ATMi 1x 2',2,9;
    'ATMi 1x 3',3,9;
    'ATMi 1x 4',4,9;
    'ATMi 1x 5',5,9;
    'ATMi 1x 6',6,9;
    %'ATMi 1x 7',7,9; %low sample
    %'ATMi 1x 8',8,9; %low sample
    %}
}
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'Units','pixels');                %sets screensize units by pixels
screendims = get(0,'ScreenSize');       %get screensize in pixels
screenx = screendims(3);
screeny = screendims(4);

condnum = size(conditions,1);
samplesize = zeros(condnum);
%namecount = cell(condnum,1);
namecount = cell(2,1);
totaltime=40;
%bmin = 6;       %S-G2 length
%bmax = 16;      %S-G2 length
bmin = 0;       %G1 length      drugs: 0   siRNA: 0
bmax = 15;      %G1 length      drugs: 20   siRNA: 10
%bmax = 26;
%bmin = 0.5;       %whole: 0.9    1st: 0.7    last: 0.8     %G1: 0.4
%bmax = 1.2;     %whole: 1.7    1st: 1.5    last: 2.0     %G1: 1.2

bstep = (bmax-bmin)/(bmax*80);
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
if stairview || scatterview
    stairfig = figure;
    %set(stairfig,'Position',[round(0.4*screenx) round(0.3*screeny) round(0.3*screenx) round(0.4*screeny)]); %1x1
    %set(stairfig,'color','w','PaperPosition',[0 0 6 4]); %cumulative stairs (1x1)
    %set(stairfig,'Position',[round(0.4*screenx) round(0.3*screeny) round(0.6*screenx) round(0.4*screeny)]); %3x1
    %set(stairfig,'color','w','PaperPosition',[0 0 14 4]); %cumulative stairs (3x1)
    
    set(stairfig,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]); %3x3
    set(stairfig,'color','w','PaperPosition',[0 0 18 12]); %cumulative stairs (3x4)
    
    %set(stairfig,'color','w','PaperPosition',[0 0 18 12]); %cumulative stairs (3x3)
    %set(stairfig,'color','w','PaperPosition',[0 0 12 8]); %cumulative stairs (3x2)
    palette = 'krgc';
end

for i = 1:condnum
    hold on;
    row = cell2mat(conditions(i,2));
    col = cell2mat(conditions(i,3));
    cd ..\Analysis
    [firstmitosis,g1length,g1lengthlineage,sg2length,angieg1,angiemean] = cdt1postdrugrepeat_pastcorrelate_cdt1only(row,col,bmax);
    
    %%% value selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    valuesfm = firstmitosis(:);
    values = g1length(:);
    %values = g1lengthlineage;
    %values = sg2length(:);
    %values = angieg1(:);
    %values = angiemean(:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    samplesize(i) = numel(values);
    cd ..\Functions
    if histview
        figure(histfig);
        %subaxis(3,1,i,'SV',0.05,'MR',0.05,'MT',0.1,'MB',0.1,'PB',0.07); %5x4
        %subaxis(4,1,i,'MR',0.05,'MT',0.1,'MB',0.1,'PB',0.05); %Cdt 4x1
        %subaxis(3,1,i,'SV',0.08,'ML',0.2,'MR',0.05); %Cdt 3x1
        %subaxis(2,1,i,'SV',0.08,'ML',0.2,'MR',0.05,'PT',0.1); %Cdt 2x1        
        hist(values,bin);
        noentrypercent = 100*noentry/(noentry+entry);
        xlabel('CDK2-Cyclin A activity');
        %xlabel('S-M time length (hrs)');
        %xlabel([sprintf('%0.2f',noentrypercent),'% arrested']);
        ylabel('cumulative distribution (%)');
    end
    if stairview
        n_elements = histc(values,bin);
        c_elements = cumsum(n_elements);
        n_elements = 100*n_elements/max(n_elements);
        c_elements = 100*c_elements/max(c_elements);
        c_elements(c_elements==max(c_elements))=99;         %stupid, but can't see top line otherwise
        if i==1
            control = c_elements;
        else
            %subaxis(1,1,i-1,'SV',0.05,'ML',0.05,'MR',0.05,'PB',0.05); %Cdt 1x1
            %subaxis(1,3,i-1,'SV',0.05,'ML',0.05,'MR',0.05,'PB',0.05); %Cdt 3x1
            %subaxis(3,2,i-1,'ML',0.05,'MR',0.05,'MT',0.01,'MB',0.03); %Cdt 3x2
            subaxis(3,3,i-1,'ML',0.05,'MR',0.05,'MT',0.01,'MB',0.03); %Cdt 3x3
            %subaxis(3,4,i-1,'ML',0.05,'MR',0.05,'MT',0.01,'MB',0.03); %Cdt 3x4
            stairs(bin,control,'k','linewidth',2); hold on;
            stairs(bin,c_elements,'r','linewidth',2);
            %stairs(bin,n_elements,palette(i),'linewidth',2);
            xlim([bmin bmax]);
            k=0;
            for j=[1 i]
                k=k+1;
                namecount{k} = [char(conditions(j,1)),' (n = ',mat2str(samplesize(j)),')'];
            end
            legend(char(namecount(:)),'location','northwest');
        end
    end
    if scatterview
        %subaxis(1,3,i,'SV',0.05,'ML',0.05,'MR',0.05,'PB',0.05); %Cdt 3x1
        %subaxis(3,2,i,'ML',0.05,'MR',0.05,'MT',0.01,'MB',0.03); %Cdt 3x2
        subaxis(3,3,i,'ML',0.05,'MR',0.05,'MT',0.03,'MB',0.03); %Cdt 3x3
        %subaxis(3,4,i,'ML',0.05,'MR',0.05,'MT',0.05,'MB',0.03); %Cdt 3x4
        %plot(g1lengthlineage(:,1),g1lengthlineage(:,2),'.','markersize',20);
        plot(valuesfm,values,'.','markersize',20);
        title(char(conditions(i,1)));
        %axis([bmin 20 bmin bmax]);
        axis([0 totaltime-bmax bmin bmax-1]);
    end
end
saveas(gcf,'h:\Downloads\Fig.jpg');
%close(histfig);
cd ..\Analysis; %return to this directory