clear; close all;
cd ..\Functions; %change directory for function calls
histview = 1;
stairview = 0;
scatterview = 0;

%%%%%%%% define conditions %%%%%%%%%%

conditions = {
    %'Media',3,1;
    'DMSO',3,2;
    %'Gefitinib 10uM',4,1;        %full block
    'Harmine 5uM',4,10;         %late block
    %{
    'Harmine 0.5uM',3,10;         %late block
    'Harmine 5uM',4,10;         %late block
    'Ponatinib 10nM',3,12;       %accelerate
    'Ponatinib 100nM',4,12;       %accelerate
    'U73122 0.5uM',3,3;           %early accelerate
    %'U73122 5uM',4,3;           %early accelerate + late block
    %'Forskolin 1uM',3,5;        %delay
    'Forskolin 10uM',4,5;        %delay
    'Thapsigargin 10nM',3,4;     %big delay
    %'Rapamycin 10nM',3,6;        %delay + early delay
    'Rapamycin 100nM',4,6;        %delay + early delay
    %'U0126 0.1uM',3,7;            %delay + late delay
    'U0126 1uM',4,7;            %delay + late delay
    %'LY294002 0.1uM',3,8;         %synchronizer
    'LY294002 1uM',4,8;         %synchronizer
    %'Ruboxistaurin',3,9;    %early delay
    'Nutlin 1uM',3,11;          %full block
    %'Nutlin 10uM',4,11;          %full block
    
    %}
    %{
    'Harmine 0.5uM',3,10;         %late block
    'Harmine 5uM',4,10;           %late block
    'U73122 0.5uM',3,3;           %early accelerate
    'LY294002 0.1uM',3,8;         %synchronizer
    'Ponatinib 10nM',3,12;       %accelerate
    'Ponatinib 100nM',4,12;       %accelerate
    %}
}
%}
%{
conditions = {
    'Non-specific siRNA',[3 4 5],1; %6.5-30.5hr
    'p21 siRNA',[3 4 5],11;         %6.5-30.5hr
    %'Non-specific siRNA',[5 6 7],1; %24-48hr
    %'p21 siRNA',[5 6 7],11;         %24-48hr
    %'Non-specific siRNA',[5 6 7],1; %48-72hr
    %'p21 siRNA',[5 6 7],10;         %48-72hr    
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
%bmin = 6;       %S-G2 length
%bmax = 16;      %S-G2 length
bmin = 0;       %G1 length      drugs: 0   siRNA: 0
bmax = 16;      %G1 length      drugs: 20   siRNA: 10
%bmin = 0.5;       %whole: 0.9    1st: 0.7    last: 0.8     %G1: 0.4
%bmax = 1.2;     %whole: 1.7    1st: 1.5    last: 2.0     %G1: 1.2

bstep = (bmax-bmin)/(bmax*1);
bin = [bmin:bstep:bmax];

%%% set up print options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set(gcf,'Position',[round(0.4*screenx) round(0.3*screeny) round(0.3*screenx) round(0.3*screeny)]); %1x1
%set(gcf,'color','w','PaperPosition',[0 0 6 4]); %cumulative stairs (1x1)
%set(gcf,'Position',[round(0.4*screenx) round(0.3*screeny) round(0.6*screenx) round(0.4*screeny)]); %3x1
%set(gcf,'color','w','PaperPosition',[0 0 14 4]); %cumulative stairs (3x1)

set(gcf,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]); %3x3
set(gcf,'color','w','PaperPosition',[0 0 18 12]); %cumulative stairs (3x4)

%set(gcf,'color','w','PaperPosition',[0 0 18 12]); %cumulative stairs (3x3)
%set(gcf,'color','w','PaperPosition',[0 0 12 8]); %cumulative stairs (3x2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:condnum
    hold on;
    row = cell2mat(conditions(i,2));
    col = cell2mat(conditions(i,3));
    cd ..\Analysis
    [firstmitosis,g1length,g1lengthlineage,sg2length,angieg1,angiemean] = cdt1postdrugrepeat_pastcorrelate(row,col,bmax);
    
    %%% value selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %values = firstmitosis(:);
    values = g1length(:);
    %values = g1lengthlineage;
    %values = sg2length(:);
    %values = angieg1(:);
    %values = angiemean(:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    samplesize(i) = numel(values);
    cd ..\Functions
    if histview || stairview
        bn_elements = histc(values,[0,6,14,bmax+1]);
        bn_elements = 100*bn_elements/sum(bn_elements);
        bn_elements
        kn_elements = ksdensity(values,bin);
        n_elements = histc(values,bin);
        c_elements = cumsum(n_elements);
        c_elements = 100*c_elements/max(c_elements);
        c_elements(c_elements==max(c_elements))=99;         %stupid, but can't see top line otherwise
        if i==1
            k_control = kn_elements;
            control = c_elements;
        else
            %subaxis(1,3,i-1,'SV',0.05,'ML',0.05,'MR',0.05,'PB',0.05); %Cdt 3x1
            %subaxis(3,2,i-1,'ML',0.05,'MR',0.05,'MT',0.01,'MB',0.03); %Cdt 3x2
            %subaxis(3,3,i-1,'ML',0.05,'MR',0.05,'MT',0.01,'MB',0.03); %Cdt 3x3
            subaxis(3,4,i-1,'ML',0.05,'MR',0.05,'MT',0.01,'MB',0.03); %Cdt 3x4
            if histview
                plot(bin,k_control,'k','linewidth',2); hold on;
                plot(bin,kn_elements,'r','linewidth',2);
            elseif stairview
                stairs(bin,control,'k','linewidth',2); hold on;
                stairs(bin,c_elements,'r','linewidth',2);
            end
            xlim([bmin bmax]);
            k=0;
            for j=[1 i]
                k=k+1;
                namecount{k} = [char(conditions(j,1)),' (n = ',mat2str(samplesize(j)),')'];
            end
            legend(char(namecount(:)),'location','northeast');
        end
    end
    if scatterview
        %subaxis(1,3,i-1,'SV',0.05,'ML',0.05,'MR',0.05,'PB',0.05); %Cdt 3x1
        %subaxis(3,2,i-1,'ML',0.05,'MR',0.05,'MT',0.01,'MB',0.03); %Cdt 3x2
        %subaxis(3,3,i-1,'ML',0.05,'MR',0.05,'MT',0.01,'MB',0.03); %Cdt 3x3
        subaxis(3,4,i,'ML',0.05,'MR',0.05,'MT',0.05,'MB',0.03); %Cdt 3x4
        plot(g1lengthlineage(:,1),g1lengthlineage(:,2),'.','markersize',20);
        title(char(conditions(i,1)));
        axis([bmin 20 bmin bmax]);
    end
end
saveas(gcf,'h:\Downloads\Fig.jpg');
%close(gcf);
cd ..\Analysis; %return to this directory