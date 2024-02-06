%function [mitosistime,g1length] = getG1length(row,col,g1max,drugspike)
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
leavein = ~tempzeros;   %ones

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
        [~,g1length(i,k)] = max(cdt1(i,curpoint+1:nextpoint));
        sphasestart(i,k) = curpoint+g1length(i,k);
        boundaryflag(i,k) = truemax(sphasestart(i,k),curpoint,nextpoint);
        k = k+1;
    end
end
%%%% gate out cells that fall outside our time limits (to avoid biasing) %%%%%%%%%%%%%%%%%%%%%%
leavein(mitosesafterspike==0)=0;
leavein(mitosesafterspike>=maxmitosis)=0;
leavein(logical(boundaryflag))=0;
g1length(g1length>g1max)=g1max;     %clip value to represent arrest
%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 1;
if panelvisualize
    %for i=1:samplesize
    for i=1:72    
        figure(ceil(i/24));             %only 24 plots per figure
        subplot(4,6,mod(i-1,24)+1);
        cdt1smooth = smooth(cdt1(i,:));
        angiesmooth = smooth(angie(i,:));
        [haxes,hline1,hline2] = plotyy(xtime,angiesmooth,xtime,cdt1smooth);
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
            Anaphases = mitosesafterspike(i,markable)+1;
            plot(Anaphases,angiesmooth(Anaphases),'ro','markerfacecolor','r');
        end
        axes(haxes(2)); hold on;
        axis([0 numframes ymincdt1 ymaxcdt1]);
        set(gca,'YTick',ymincdt1:ystepcdt1:ymaxcdt1);
        set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
        if markable
            SphaseEntries = sphasestart(i,markable);
            plot(SphaseEntries,cdt1smooth(SphaseEntries),'ko','markerfacecolor','k');
        end
    end
end
%%%%%%%%%%%%%%%%%% Remove bad data and sort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mitosistime = mitosesafterspike(logical(leavein));
g1length = g1length(logical(leavein));       %goodg1length = g1length.*countg1;

cd ..\Analysis; %return to this directory