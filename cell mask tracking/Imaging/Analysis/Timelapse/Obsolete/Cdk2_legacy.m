clear; close all;
cd('H:\Documents\Timelapse\Code\Development\Functions\'); %change directory for function calls
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130715\';
stairview = 1;
timeview = 0;
timepercentview = 0;
correlationview = 0;

%%%%%%%% define conditions %%%%%%%%%%
conditions = {
    'DMSO All',1:8,10;
}
%{
conditions = {
    'DMSO All',1:8,10;
    'Media All',1:8,11;
    'MK-2206 All',[1 2 3 4 5 6 8],4;
    'SB203580 All',1:8,6;
    'Harmine All',1:8,8;
    'ATMi All',[2 3 4 5],9;
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

numframes = 231;
drugspike = 0;
framesperhr = 5;
measurelimit = 15;    %this is the longest measurement to consider
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

%%% figure settings %%%
stairfig = figure;
set(stairfig,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]); %3x3
set(stairfig,'color','w','PaperPosition',[0 0 12 8]); %cumulative stairs (2x3)
palette = 'krgc';
%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:condnum
    hold on;
    row = cell2mat(conditions(i,2));
    col = cell2mat(conditions(i,3));
    [goodsenescence,badsenescence] = tracestats_senescence(row,col,site,path);
    [mitosis,IMT] = tracefeatures_IMT(row,col,measurelimit,drugspike);
    %[mitosis,g1entry,g0toolong] = tracefeatures_cdk2(row,col,measurelimit,drugspike);
    %%% value selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    times = mitosis/framesperhr;
    %times = (mitosis+g1entry)/framesperhr;
    %times = g1entry/framesperhr;
    
    %values = sentry/framesperhr;     %g1 length
    %values = (sentry-g1entry)/framesperhr;      %g1 rise length
    %values = g1entry/framesperhr;      %g1 time till rise ('latency')
    values = g0entry/framesperhr;      %g1 fall length
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
    quiescentcells = g1toolong;   %set to g1toolong or g0toolong
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
cd('H:\Documents\Timelapse\Code\Development\Analysis\Timelapse\'); %return to this directory