function [outputdata] = cdt1postdrug(row,col)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timescape\20120807_12drugs\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120202_6p5-30p5hr\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120209_24-48hr\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120115_48-72hr\';
%%%%%%%%%%%% Experiment Settings %%%%%%%%
nocall = 0;
if nocall
    row = [4];
    col = [10];
end
site = [0 1];
framesperhr = 5;
drugspike = 90;          %set to zero if no drug
yminangie = 0.4; ymaxangie = 2; ystepangie = 0.2;
ymincdt1 = 0; ymaxcdt1 = 3; ystepcdt1 = 0.5;
signalwidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[signal1,signal2,mitoses] = combinedata(row,col,site,path,drugspike);
[signal1,signal2,mitoses] = combinedata_realbirths(row,col,site,path);
tracelist = 1:size(signal2,1);
%%%%%%% remove cells that never mitose %%%%%%%%%%%%%%%%
%fprintf(['count after removing duplicates: ',num2str(length(tracelist)),'\n']);
tracelist(ismember(tracelist,find(max(mitoses,[],2)==0))) = [];
samplesize = length(tracelist);
angie = signal1(tracelist,:);
cdt1 = signal2(tracelist,:);
xtime = 1:size(angie,2);
numframes = size(cdt1,2);
cdt1time = zeros(samplesize,1);
cdt1timeactual = zeros(samplesize,1);
goodmitoses = mitoses(tracelist,:);
%%%% mark mitoses and find first mitosis after drug spike %%%%%
mitosisafterspike = zeros(samplesize,1);    %1st mitosis after drugspike
nextmitosis = zeros(samplesize,1);          %2nd mitosis after drugspike
leaveout = zeros(samplesize,1);
boundaryflag = zeros(samplesize,1);         %flag for bad cdt1 maximum (sides don't fall off)
sg2length = zeros(samplesize,1);            %time from cdt1 degradation to next mitosis
ccend = zeros(samplesize,1);                %end of cell cycle of interest
for i=1:samplesize
    mitosistimes = goodmitoses(i,:);
    mitosistimes = mitosistimes(mitosistimes>0);
    mitosistimes = sort(mitosistimes);                  %get time window between drugspike & next mitosis
    temp = find(mitosistimes>=drugspike,1);
    if temp
        mitosisafterspike(i) = mitosistimes(temp);          %record mitosis
        if temp < length(mitosistimes)                  %if there is a subsequent mitosis (only relevant for pre-spike mitoses)
            nextmitosis(i) = mitosistimes(temp+1);
        else
            nextmitosis(i) = 0;
        end
    else
        mitosisafterspike(i) = 0;
    end
end

for i=1:samplesize
    if mitosisafterspike(i)==0
        continue        %skip to next iteration
    end
    [~,cdt1timepre] = max(cdt1(i,drugspike:mitosisafterspike(i)));
    actualtimepre = drugspike+cdt1timepre-1;
    boundaryflag(i) = truemax(cdt1(i,:),actualtimepre,numframes,drugspike);
    if ~boundaryflag(i)    %true max btwn drugspike & first mitosis
        cdt1time(i) = cdt1timepre;
        cdt1timeactual(i) = actualtimepre;
        sg2length(i) = mitosisafterspike(i)-cdt1timeactual(i);
        ccend(i) = mitosisafterspike(i);
    else
        if nextmitosis(i)
            [~,cdt1time(i)] = max(cdt1(i,mitosisafterspike(i):nextmitosis(i)));
            cdt1timeactual(i) = mitosisafterspike(i)+cdt1time(i)-1;
            boundaryflag(i) = truemax(cdt1(i,:),cdt1timeactual(i),numframes,mitosisafterspike(i));
            sg2length(i) = nextmitosis(i)-cdt1timeactual(i);
            ccend(i) = nextmitosis(i);
        end
    end
end
%%%% gate out cells that fall outside our time limits (to avoid biasing) %%%%%
leaveout(ccend==0) = 1;         %ccend=0 if no valid cell cycle end
leaveout(boundaryflag==1) = 1;  %boundaryflag=1 if max cdt1 doesn't drop off on both sides
%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 0;
if panelvisualize
    for i=1:samplesize
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
        markable = ccend(i)~=0 && boundaryflag(i)==0;
        if markable
            plot(ccend(i),angiesmooth(ccend(i)),'ro','markerfacecolor','r');
        end
        axes(haxes(2)); hold on;
        axis([0 numframes ymincdt1 ymaxcdt1]);
        set(gca,'YTick',ymincdt1:ystepcdt1:ymaxcdt1);
        set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
        if markable
            plot(cdt1timeactual(i),cdt1smooth(cdt1timeactual(i)),'ko','markerfacecolor','k');
        end
    end
end
%%% Remove bad data and sort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodindex = leaveout==0;
sg2length = sg2length/framesperhr;
outputdata = sg2length(goodindex);
%%% View Summary Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summarygraph = 0;
if summarygraph
    plot(outputdelay,outputdata,'.','markersize',20); %turn on to visualize
    xlabel('Time from Drug Spike --> Mitosis (Hrs)','Fontsize',20);
    ylabel('Time from mitosis --> Cdt1 destruction (Hrs)','Fontsize',20);
    saveas(gcf,'h:\Downloads\Fig.jpg');
end
cd ..\Analysis; %return to this directory
end

function [boundarycheck] = truemax(cdt1,actualtime,numframes,leftedge)
lowertimebound = actualtime-20;
uppertimebound = actualtime+20;
%%% check left side %%%%%%%%%%%%%%%%%
dataedge = lowertimebound<0;
stillrising = actualtime==leftedge && cdt1(max([leftedge-3,1]))>cdt1(leftedge);
slowfall = cdt1(actualtime)-cdt1(lowertimebound)<0.3;
traceedge = cdt1(actualtime-5)==-10000;
lowercheck = dataedge || stillrising || slowfall || traceedge;
%%% check right side %%%%%%%%%%%%%%%%
dataedge = uppertimebound>numframes;
slowfall = cdt1(actualtime)-cdt1(uppertimebound)<0.3;
traceedge = cdt1(actualtime+5)==-10000;
uppercheck = dataedge || slowfall || traceedge;

boundarycheck = lowercheck || uppercheck;
end
