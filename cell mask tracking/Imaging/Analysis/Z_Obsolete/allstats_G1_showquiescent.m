clear; close all;
cd ..\Functions; %change directory for function calls
stairview = 0;
timeview = 1;
timepercentview = 0;
correlationview = 0;

%%%%%%%% define conditions %%%%%%%%%%

conditions = {
    'DMSO All',1:8,10;
    'Media All',1:8,11;
    'MK-2206 All',[1 2 3 4 5 6 8],4;
    'SB203580 All',1:8,6;
    'Harmine All',1:8,8;
    'ATMi All',[2 3 4 5],9;
}
%}
%{
conditions = {
    'DMSO All',1:6,10;
    %{
    'DMSO 1',1,10;
    %'DMSO 2',2,10;
    'DMSO 3',3,10;
    'DMSO 4',4,10;
    %'DMSO 5',5,10;
    %'DMSO 6',6,10;
    %'DMSO 7',7,10;
    %'DMSO 8',8,10;
    %}
    
    %{
    'MK-2206 1uM 1',1,4;
    'MK-2206 1uM 2',2,4;
    'MK-2206 1uM 3',3,4;
    'MK-2206 1uM 4',4,4;
    'MK-2206 1uM 5',5,4;
    'MK-2206 1uM 6',6,4;
    %'MK-2206 1uM 7',7,4; %data missing
    'MK-2206 1uM 8',8,4;
    %}
    %{
    'SB203580 10uM 1',1,8;
    'SB203580 10uM 2',2,8;
    'SB203580 10uM 3',3,8;
    'SB203580 10uM 4',4,8;
    'SB203580 10uM 5',5,8;
    'SB203580 10uM 6',6,8;
    'SB203580 10uM 7',7,8;
    'SB203580 10uM 8',8,8;
    %}
    %{
    'Harmine 1uM  1',1,8;
    'Harmine 1uM  2',2,8;
    'Harmine 1uM  3',3,8; %low sample
    'Harmine 1uM  4',4,8;
    'Harmine 1uM  5',5,8;
    'Harmine 1uM  6',6,8; %low sample
    'Harmine 1uM  7',7,8; %low sample
    'Harmine 1uM 8',8,8; %low sample
    %}
    
    %{
    %'ATMi 1x 1',1,9;
    'ATMi 1x 2',2,9;
    'ATMi 1x 3',3,9;
    'ATMi 1x 4',4,9;
    'ATMi 1x 5',5,9;
    %'ATMi 1x 6',6,9;
    %'ATMi 1x 7',7,9; %low sample
    %'ATMi 1x 8',8,9; %low sample
    %}
}
%}
%{
conditions = {
    %'Media',3,1;
    'DMSO',3,2;
    %'Gefitinib 10uM',4,1;        %full block
    
    'Harmine 0.5uM',3,10;         %late block
    %'Harmine 5uM',4,10;         %late block
    'Ponatinib 10nM',3,12;       %accelerate
    %'Ponatinib 100nM',4,12;       %accelerate
    'U73122 0.5uM',3,3;           %early accelerate
    %'U73122 5uM',4,3;           %early accelerate + late block
    'Forskolin 1uM',3,5;        %delay
    'Forskolin 10uM',4,5;        %delay
    'Thapsigargin 10nM',3,4;     %big delay
    'Rapamycin 10nM',3,6;        %delay + early delay
    'Rapamycin 100nM',4,6;        %delay + early delay
    
    'U0126 0.1uM',3,7;            %delay + late delay
    'U0126 1uM',4,7;            %delay + late delay
    %'LY294002 0.1uM',3,8;         %synchronizer
    %'LY294002 1uM',4,8;         %synchronizer
    'Ruboxistaurin',3,9;    %early delay
    'Nutlin 1uM',3,11;          %full block
    %'Nutlin 10uM',4,11;          %full block
    
    
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
samplesize = zeros(condnum,1);
qper = zeros(condnum,1);
namecount = cell(2,1);

numframes = 240;
drugspike = 40;
framesperhr = 5;
measurelimit = 20;    %this is the longest measurement to consider
measuremargin = 10;   %this is the margin needed to identify cdt1 peak or hdhb rise
longthresh = 4;    %length of latency considered 'quiescent' in timepercentview

bstep = 1/framesperhr;
bmin = 0;
bmax = measurelimit;
%bmax = 15;
bin = [bmin:bstep:bmax];
%tmin = drugspike/framesperhr;
%tmax = (numframes-g1max*framesperhr)/framesperhr;
tmin = 0;
tmax = (numframes-drugspike-measurelimit*framesperhr-measuremargin)/framesperhr;

if stairview || timeview || correlationview || timepercentview
    stairfig = figure;
    set(stairfig,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]); %3x3
    set(stairfig,'color','w','PaperPosition',[0 0 12 8]); %cumulative stairs (2x3)
    palette = 'krgc';
end

for i = 1:condnum
    hold on;
    row = cell2mat(conditions(i,2));
    col = cell2mat(conditions(i,3));
    cd ..\Analysis
    [mitosis,g0entry,g1entry,sentry,g1toolong,g0toolong] = getG1features(row,col,measurelimit*framesperhr,drugspike);
    
    %%% value selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    times = mitosis/framesperhr;
    %times = (mitosis+g1entry)/framesperhr;
    %times = (mitosis+g0entry)/framesperhr;
    
    %values = sentry/framesperhr;     %g1 length
    %values = (sentry-g1entry)/framesperhr;      %g1 rise length
    values = g1entry/framesperhr;      %g1 time till rise ('latency')
    %values = g0entry/framesperhr;      %g1 fall length
    %values = (g1entry-g0entry)/framesperhr;    %g0 trough length
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% time screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    timescreen = 0;
    if timescreen
        timepoint = 7;  %number of hours after drugspike before which data is to be ignored
        timepoint = timepoint + drugspike/framesperhr;
        tooearly = find(times<timepoint);
        times(tooearly)=[];
        values(tooearly)=[];
        values2(tooearly)=[];
    end
    
    %%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    cd ..\Functions
    quiescentcells = g0toolong;   %set to g1toolong or g0toolong
    if stairview
        samplesize(i) = numel(values);
        quiescentcount = sum(sum(quiescentcells));
        values(quiescentcells==1)=[];
        %qper(i) = round(100*quiescentcount/samplesize(i));
        n_elements = histc(values,bin);
        g1risemeasured = 1;     %set to 1 if measuring G1 rise length
        if g1risemeasured
            samplesize(i)=numel(values);    %re-calc to only count valid traces
        else
            n_elements(end) = n_elements(end)+quiescentcount;     %so we can represent the quiescent cells
        end
        c_elements = cumsum(n_elements);
        %n_elements = 100*n_elements/max(n_elements);
        c_elements = 100*c_elements/max(c_elements);
        c_elements(c_elements==max(c_elements))=99;         %stupid, but can't see top line otherwise
        if i==1
            control = c_elements;
        else
            subaxis(2,3,i-1,'ML',0.1,'MR',0.03,'MT',0.02,'MB',0.1);
            stairs(bin,control,'k','linewidth',2); hold on;
            stairs(bin,c_elements,'r','linewidth',2);
            %stairs(bin,n_elements,palette(i),'linewidth',2);
            xlim([bmin bmax]);
            k=0;
            for j=[1 i]
                k=k+1;
                namecount{k} = [char(conditions(j,1)),' (n=',mat2str(samplesize(j)),')'];
            end
            legend(char(namecount(:)),'location','southeast');
        end
    end
    
    if timeview
        times = times-drugspike/framesperhr;
        times(quiescentcells==1)=[];
        values(quiescentcells==1)=[];
        [times,idx]=sort(times);
        values=values(idx);
        tbin=floor(times);
        
        subaxis(2,3,i,'ML',0.1,'MR',0.05,'MT',0.05,'MB',0.1,'SV',0.1);
        plot(times,values,'.','markersize',8);   %size 12 for presentations, 8 for screen
        title(char(conditions(i,1)));
        timemin=tmin; timemax=tmax;
        %timemin=floor(min(times)); timemax=ceil(max(times));
        axis([timemin timemax bmin bmax]);
        hold on;   
        
        continuoustrend=0; %set to 1 for spline fitting, otherwise 0 for bar averaging
        if continuoustrend
            discretetimes=max(tbin)*framesperhr;
            filtcoef = round(5*length(times)/discretetimes);   %number of values per each timestep + 2 prior and 2 post
            fprintf('%0.0f\n',filtcoef);
            filtfiltvec = (1/filtcoef)*ones(1,filtcoef);
            times = [-1*times(end:-1:1);times;times+times(end)];
            values = [values(end:-1:1);values;values(end:-1:1)];
            trend=filtfilt(filtfiltvec,1,values);
            line(times,trend,'Color','r','Linewidth',3);            
        else
            avgvals=zeros(size(values));
            tbinunique=unique(tbin);
            tbinnum=length(tbinunique);
            for cc=1:tbinnum
                vc=tbinunique(cc);
                tcidx=find(tbin==vc);
                tcavg=mean(values(tcidx));
                avgvals(tcidx)=tcavg;
            end
            scatter(times,avgvals,20,'ro','markerfacecolor','r');  
        end
    end
    
    if timepercentview
        %quiescentcells(values>longthresh)=1;  %de-comment to classify quiescence by a shorter time
        times = times-drugspike/framesperhr;
        [times,idx]=sort(times);
        values=values(idx);
        quiescentcells=quiescentcells(idx);
        values(quiescentcells==1)=0;  %avoid counting latency for traces with g1length>g1max
        tbin=floor(times); %bin by hour
        %tbin=times;         %no bin
        tbinunique=unique(tbin);
        tbinnum=length(tbinunique);
        binqper=zeros(tbinnum,1);
        for cc=1:tbinnum
            vc=tbinunique(cc);
            tcidx=find(tbin==vc);
            numtotal=length(tcidx);
            numquiescent=sum(quiescentcells(tcidx));
            binqper(cc)=round(100*numquiescent/numtotal);  %returns percent quiescent
            %numlong=sum(values(tcidx)>=longthresh);
            %binqper(cc)=round(100*numlong/(numtotal-numquiescent));
        end
        subaxis(2,3,i,'ML',0.1,'MR',0.05,'MT',0.05,'MB',0.1,'SV',0.1);
        bar(tbinunique,binqper,0.4,'g','EdgeColor','k');
        
        %{
        %%% apply height shading %%%
        h=bar(tbinunique,binqper);
        ch=get(h,'Children');
        fvd=get(ch,'Faces');
        fvcd=get(ch,'FaceVertexCData');
        [zs,izs]=sortrows(binqper,1);
        k=128;
        colormap(summer(k));
        shading interp
        for ix=1:tbinnum
            color=floor(k*zs(ix)/tbinnum);
            row=izs(ix);
            fvcd(fvd(row,1))=1;
            fvcd(fvd(row,4))=1;
            fvcd(fvd(row,2))=color;
            fvcd(fvd(row,3))=color;
        end
        set(ch,'FaceVertexCData',fvcd);
        set(ch,'EdgeColor','k');
        %}
        
        title(char(conditions(i,1)));
        axis([-0.5 tmax bmin 100]);
    end
    
    if correlationview
        b1min = 0; b1max = 15;
        b2min = 0; b2max = 15;
        subaxis(2,3,i,'ML',0.1,'MR',0.05,'MT',0.05,'MB',0.1,'SV',0.1);
        plot(values2,values,'.','markersize',8);
        %dscatter(values2,values,'PLOTTYPE','contour');
        title(char(conditions(i,1)));
        axis([b2min b2max b1min b1max]);
    end
end
saveas(gcf,'h:\Downloads\Fig.jpg');
%close(histfig);
cd ..\Analysis; %return to this directory