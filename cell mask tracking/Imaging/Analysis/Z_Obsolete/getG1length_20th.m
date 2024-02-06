%function [mitosistime,g1length,g1droplength,g1botlength,g1riselength] = getG1length_test(row,col,g1max,drugspike)
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
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120202_6p5-30p5hr\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120209_24-48hr\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120115_48-72hr\';
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
%%%%%%%%%%%% Experiment Settings %%%%%%%%%%%%%%%%
nocall = 0;
if nocall
    row = [];
    col = [];
end
site = [1];
%%%%%%%%%%%%%% viewer settings %%%%%%%%%%%%%%%%%%
yminangie = 0.3; ymaxangie = 2.5; ystepangie = 0.2;
ymincdt1 = 0; ymaxcdt1 = 600; ystepcdt1 = 50;
signalwidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[signal1,signal2,mitoses] = combinedata_realbirths(row,col,site,path);
tracelist = 1:size(signal1,1);  %was signal2
tracelist(ismember(tracelist,find(max(mitoses,[],2)==0))) = [];  %remove any traces with no mitoses
angie = signal1(tracelist,:);
cdt1 = signal2(tracelist,:);
angiesmooth = zeros(size(angie));
cdt1smooth = zeros(size(cdt1));
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
threshpercent = 0.2;
g1botstartrel = tempzeros;
g1botstartactual = tempzeros;
g1botendrel = tempzeros;
g1botendactual = tempzeros;
leavein = ~tempzeros;   %ones

for i=1:samplesize
    temp = goodmitoses(i,:);
    temp = sort(temp(temp>drugspike));
    mitosesafterspike(i,1:length(temp)) = temp;
end


for i=1:samplesize
    angiesmooth(i,:) = smooth(angie(i,:));
    cdt1smooth(i,:) = smooth(cdt1(i,:));
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
        %[~,g1length(i,k)] = max(smooth(cdt1(i,curpoint+1:nextpoint)));
        [~,g1length(i,k)] = max(cdt1smooth(i,curpoint+1:nextpoint));
        sphasestart(i,k) = curpoint+g1length(i,k);
        boundaryflag(i,k) = truemax(sphasestart(i,k),curpoint,nextpoint);   %check if likely error
        if boundaryflag(i,k)==0
            
            %g1stretch = smooth(angie(i,curpoint+1:curpoint+g1length(i,k)));
            g1stretch = angiesmooth(i,curpoint+1:curpoint+g1length(i,k));
            
            %{
            tempsort = sort(g1stretch);
            g1lowliers = find(g1stretch<tempsort(2)-0.05);
            if g1lowliers
                g1stretchprev = [g1stretch(2);g1stretch(1:end-1)];
                g1stretchpost = [g1stretch(2:end);g1stretch(end-1)];
                g1stretchinterpolate = (g1stretchprev+g1stretchpost)/2;
                g1stretch(g1lowliers) = g1stretchinterpolate(g1lowliers);
            end
            %}
            
            lowg1thresh = min(g1stretch)+range(g1stretch)*threshpercent;
            lowg1 = find(g1stretch<=lowg1thresh);
            
            %{
            g1stretchdiff = diff(g1stretch); g1stretchdiff(end)=0; g1stretchdiff=[g1stretchdiff;0]; %append 0
            g1upslope = g1stretchdiff>-0.02;
            g1downslope = g1stretchdiff<0.02;
            g1botstartrel(i,k) = find(lowg1 & g1upslope,1);
            g1botendrel(i,k) = find(lowg1 & g1downslope,1,'last')+1;
            %}
            
            g1botstartrel(i,k) = lowg1(1);
            g1botendrel(i,k) = lowg1(end);
            
            g1botstartactual(i,k) = curpoint+g1botstartrel(i,k);
            g1botendactual(i,k) = curpoint+g1botendrel(i,k);
        end
        k = k+1;
    end
end
g1length(g1length>g1max)=g1max;     %clip value to represent arrest
g1droplength = g1botstartrel;
g1botlength = g1botendrel-g1botstartrel;    %length of g1 trough
g1riselength = g1length-g1botendrel;        %length of g1 rise
%%%% gate out cells that fall outside our time limits (to avoid biasing) %%%%%%%%%%%%%%%%%%%%%%
leavein(mitosesafterspike==0)=0;
leavein(mitosesafterspike>=maxmitosis)=0;
leavein(logical(boundaryflag))=0;
%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 1;
if panelvisualize
    for i=1:24
        figure(ceil(i/6));             %only 24 plots per figure
        subplot(2,3,mod(i-1,6)+1);
        [haxes,hline1,hline2] = plotyy(xtime,angiesmooth(i,:),xtime,cdt1smooth(i,:));
        axes(haxes(1)); hold on;
        axis([0 numframes yminangie ymaxangie]);
        set(gca,'YTick',yminangie:ystepangie:ymaxangie);
        set(hline1,'color','b','linewidth',signalwidth);
        title(['Cell ',num2str(i)]);
        if drugspike
            line([drugspike drugspike],[yminangie ymaxangie],'Color','r');  %add vertical line when drug is added
        end
        markable = find(leavein(i,:));
        if markable
            Anaphases = mitosesafterspike(i,markable);
            plot(Anaphases,angiesmooth(i,Anaphases),'ro','markerfacecolor','r');
            G1botstarts = g1botstartactual(i,markable);
            plot(G1botstarts,angiesmooth(i,G1botstarts),'go','markerfacecolor','g');
            G1botends = g1botendactual(i,markable);
            plot(G1botends,angiesmooth(i,G1botends),'bo','markerfacecolor','b');
        end
        axes(haxes(2)); hold on;
        axis([0 numframes ymincdt1 ymaxcdt1]);
        set(gca,'YTick',ymincdt1:ystepcdt1:ymaxcdt1);
        set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
        if markable
            SphaseEntries = sphasestart(i,markable);
            plot(SphaseEntries,cdt1smooth(i,SphaseEntries),'ko','markerfacecolor','k');
        end
    end
end
%%%%%%%%%%%%%%%%%% Remove bad data and sort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mitosistime = mitosesafterspike(logical(leavein));
g1length = g1length(logical(leavein));       %goodg1length = g1length.*countg1;
g1droplength = g1droplength(logical(leavein));
g1botlength = g1botlength(logical(leavein));
g1riselength = g1riselength(logical(leavein));

cd ..\Analysis; %return to this directory