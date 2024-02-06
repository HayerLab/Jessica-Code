%function [mitosistime,g1length,g1riselength] = getG1length_rise(row,col,g1max,drugspike)
row=1;col=11;g1max=100;drugspike=40;
%%%% Plot signal after post-drugspike mitosis vs drug-mitosis delay %%%%%%%
% Author: Mingyu Chung
% Last Revision: 10/20/2012
%
% Summary: This code calculates signal values after the first mitosis
% following drugspike, and plots it over the time delay between the
% drugspike and subsequent mitosis.
%
% Use:
% Before running, execute the following code:
% 1. calcnuclei_saveimages: save cropped & processed images, and store
% features for each cell in each image
% 2. tracenuclei: track cells
% 3. getcellhistory: retrieves all desired cells.  For each cell:
%    - relevant signal (e.g. angie sensor, Ctd1, etc.)
%    - mitoses (frame for every mitosis for every cell)
% Set the 'timeafter' (in frames) for the time window to average the signal
% following post-drugspike mitosis.
% 
% This code requires:
% - combinedata(row,col,site,path,drugspike)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Directory Settings %%%%%%%%%%%%%%%%%
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timescape\20120807_12drugs\';
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
%%%%%%%%%%%% Experiment Settings %%%%%%%%%%%%%%%%
nocall = 0;
if nocall
    row = [];
    col = [];
end
site = [1];
%%%%%%%%%%%%%% viewer settings %%%%%%%%%%%%%%%%%%
y1min = 0; y1max = 2.5; y1step = 0.5;
y2min = 0; y2max = 0.6; y2step = 0.06;
signalwidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[signal1,signal2,mitoses] = combinedata_realbirths(row,col,site,path);
tracelist = 1:size(signal1,1);  %was signal2
tracelist(ismember(tracelist,find(max(mitoses,[],2)==0))) = [];  %remove any traces with no mitoses
angie = signal1(tracelist,:);
cdt1 = signal2(tracelist,:);
samplesize = length(tracelist);
numframes = size(signal1,2);
xtime = 1:numframes;
maxmitosis = numframes-g1max;  %latest possible mitosis
goodmitoses = mitoses(tracelist,:);
bufsize = 10; tempzeros = zeros(samplesize,bufsize);
mitosesafterspike = tempzeros;    %1st mitosis after drugspike
boundaryflag = tempzeros;
g1length = tempzeros;
sphasestart = tempzeros;
threshpercent = 0.15;
g1risepoint = cell(samplesize,bufsize);
leavein = ~tempzeros;   %ones
maxprev = 20;


for i=1:samplesize
    temp = goodmitoses(i,:);
    temp = sort(temp(temp>drugspike));
    mitosesafterspike(i,1:length(temp)) = temp;
end

for i=1:samplesize
    if mitosesafterspike(i,1)==0
        continue        %skip to next iteration.  any traces without postspike mitosis ignored
    end
    k = 1;
    while mitosesafterspike(i,k)
        curpoint = mitosesafterspike(i,k);
        if mitosesafterspike(i,k+1)
            nextpoint = mitosesafterspike(i,k+1);
        else
            nextpoint = numframes;
        end
        [~,g1length(i,k)] = max(smooth(cdt1(i,curpoint+1:nextpoint)));
        sphasestart(i,k) = curpoint+g1length(i,k);
        boundaryflag(i,k) = truemax(sphasestart(i,k),curpoint,nextpoint);   %check if likely error
        if boundaryflag(i,k)==0
            smoothsignal = smooth(angie(i,:));
            g1stretch = smooth(angie(i,curpoint+1:curpoint+g1length(i,k)));
            
            lowg1thresh = min(g1stretch)+range(g1stretch)*threshpercent;
            lowg1 = g1stretch<=lowg1thresh;

            %[lowg1_la,maxlabel] = bwlabel(lowg1);
            %relpos = find(lowg1_la==maxlabel);
            relpos = find(lowg1,1,'first'):find(lowg1,1,'last');
            if maxprev<length(relpos)
                relpos = relpos(end-maxprev+1:end);
            end
            
            g1slope = getslope(smoothsignal',3);
            g1curve = getcurvature(smoothsignal',6);
            g1depth = 0.5-smoothsignal';
            risefilter = smooth(g1slope*4+g1curve*2+g1depth*0.5);
            
            risevals = risefilter(relpos+curpoint);
            g1risepoint{i,k} = find(risefilter==max(risevals));
            %g1risepoint{i,k} = find(lowg1)+curpoint;
            
        end
        k = k+1;
    end
end
g1length(g1length>g1max)=g1max;     %clip value to represent arrest
%g1riselength = g1length-g1botendrel;        %length of g1 rise
%%%% gate out cells that fall outside our time limits (to avoid biasing) %%%%%%%%%%%%%%%%%%%%%%
leavein(mitosesafterspike==0)=0;
leavein(mitosesafterspike>=maxmitosis)=0;
leavein(logical(boundaryflag))=0;
%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 1;
if panelvisualize
    for i=1:20
        figure(ceil(i/20));             %only 24 plots per figure
        subaxis(4,5,mod(i-1,20)+1,'ML',0.02,'MR',0.02,'MT',0.03,'MB',0.03,'SH',0.03); %5x4
        %signal1trace = angie(i,:);
        %signal2trace = smooth(cdt1(i,:));
        %signal2trace = 0;
        
        signal1trace = smooth(angie(i,:));
        grad3 = getslope(signal1trace',3);
        acc3 = getcurvature(signal1trace',6);
        depth = 0.5-signal1trace';
        signal2trace = smooth(grad3*4+acc3*2+depth*0.5);
        
        [haxes,hline1,hline2] = plotyy(xtime,signal1trace,xtime,signal2trace);
        axes(haxes(1)); hold on;
        axis([0 numframes y1min y1max]);
        set(gca,'YTick',y1min:y1step:y1max);
        set(hline1,'color','b','linewidth',signalwidth);
        title(['Cell ',num2str(i)]);
        
        if drugspike
            line([drugspike drugspike],[y1min y1max],'Color','r');  %add vertical line when drug is added
        end
        markable = find(leavein(i,:));
        if markable
            Anaphases = mitosesafterspike(i,markable);
            plot(Anaphases,signal1trace(Anaphases),'ro','markerfacecolor','r');
            for j=1:length(markable)
                jm=markable(j);
                scatter(g1risepoint{i,jm},signal1trace(g1risepoint{i,jm}),12,'ro','markerfacecolor','r');
            end
        end

        axes(haxes(2)); hold on;
        axis([0 numframes y2min y2max]);
        set(gca,'YTick',y2min:y2step:y2max);
        set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
        line([xtime(1) xtime(end)],[0 0],'Color','k');
    end
end
%%%%%%%%%%%%%%%%%% Remove bad data and sort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mitosistime = mitosesafterspike(logical(leavein));
g1length = g1length(logical(leavein));       %goodg1length = g1length.*countg1;
%g1riselength = g1riselength(logical(leavein));

cd ..\Analysis; %return to this directory