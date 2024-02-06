function [goodfirstmitosis,goodg1length,arrested,goodsg2length,goodangieg1,goodangiemean] = cdt1postdrugrepeat(row,col)
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
%path = 'h:\Documents\Timescape\20120807_12drugs\';
path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120202_6p5-30p5hr\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120209_24-48hr\';
%path = 'h:\Documents\Timescape\Sabrina_p21_siRNA_MCF10A\20120115_48-72hr\';
%%%%%%%%%%%% Experiment Settings %%%%%%%%
nocall = 0;
if nocall
    row = [3];
    col = [2];
end
site = [0 1];
framesperhr = 5;
drugspike = 0;          %set to zero if no drug
maxfirstmitosis = 16*5;
g1window = [9 10];
%%% viewer settings %%%%%%%%%%%%%%%%%%%%%
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
g1max = numframes-drugspike - maxfirstmitosis;
cdt1time = zeros(samplesize,10);
cdt1timeactual = zeros(samplesize,10);
goodmitoses = mitoses(tracelist,:);
%%%% mark mitoses and find first mitosis after drug spike %%%%%
mitosesafterspike = zeros(samplesize,10);    %1st mitosis after drugspike
leavein = ones(samplesize,10);
countg1 = zeros(samplesize,1);
arrested = zeros(samplesize,1);
boundaryflag = zeros(samplesize,10);         %flag for bad cdt1 maximum (sides don't fall off)
boundaryflagnext = zeros(samplesize,10);
angieprevmean = zeros(samplesize,10);
angieg1 = zeros(samplesize,10);
sg2length = zeros(samplesize,10);            %time from cdt1 degradation to next mitosis
g1length = zeros(samplesize,10);
for i=1:samplesize
    temp = goodmitoses(i,:);
    temp = sort(temp(temp>drugspike));
    mitosesafterspike(i,1:length(temp)) = temp;
end

for i=1:samplesize
    if mitosesafterspike(i,1)==0
        continue        %skip to next iteration
    end
    startpoint = drugspike;
    k = 1;
    while mitosesafterspike(i,k)
        [~,cdt1time(i,k)] = max(cdt1(i,startpoint+1:mitosesafterspike(i,k)));
        cdt1timeactual(i,k) = startpoint+cdt1time(i,k);
        sg2length(i,k) = mitosesafterspike(i,k)-cdt1timeactual(i,k);
        if mitosesafterspike(i,k+1)
            nextpoint = mitosesafterspike(i,k+1);
        else
            nextpoint = numframes;
        end
        [~,g1length(i,k)] = max(cdt1(i,mitosesafterspike(i,k)+1:nextpoint));
        boundaryflagnext(i,k) = truemax(cdt1(i,:),g1length(i,k)+mitosesafterspike(i,k),numframes,mitosesafterspike(i,k));
        if mitosesafterspike(i,k)+1+g1window(2)<=numframes
            angieg1(i,k) = mean(angie(i,mitosesafterspike(i,k)+1+g1window(1):mitosesafterspike(i,k)+1+g1window(2)));
        end
        fraction = round((mitosesafterspike(i,k)-cdt1timeactual(i,k))/3);
            angieprevmean(i,k) = mean(angie(i,cdt1timeactual(i,k):mitosesafterspike(i,k)));                       %whole
            %angieprevmean(i,k) = mean(angie(i,cdt1timeactual(i,k):cdt1timeactual(i,k)+fraction));                         %1st fraction
            %angieprevmean(i,k) = mean(angie(i,cdt1timeactual(i,k)+fraction:mitosesafterspike(i,k)-fraction));      %middle fraction
            %angieprevmean(i,k) = mean(angie(i,mitosesafterspike(i,k)-fraction:mitosesafterspike(i,k)));     %last fraction
        boundaryflag(i,k) = truemax(cdt1(i,:),cdt1timeactual(i,k),numframes,startpoint);
        startpoint = mitosesafterspike(i,k);
        k = k+1;
    end
end
%%%% gate out cells that fall outside our time limits (to avoid biasing) %%%%%
leavein(mitosesafterspike==0) = 0;         %ccend=0 if no valid cell cycle end
leavein(logical(boundaryflag)) = 0;  %boundaryflag=1 if max cdt1 doesn't drop off on both sides
firstmitosis = mitosesafterspike(:,1)-drugspike;
countg1(firstmitosis>0 & firstmitosis<maxfirstmitosis) = 1;
arrested(countg1 & (boundaryflagnext(:,1) | g1length(:,1)>g1max)) = 1;
countg1(logical(arrested)) = 0;
arrested = arrested(arrested>0);
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
        markable = find(leavein(i,:));
        if markable
            plot(mitosesafterspike(i,markable),angiesmooth(mitosesafterspike(i,markable)),'ro','markerfacecolor','r');
        end
        axes(haxes(2)); hold on;
        axis([0 numframes ymincdt1 ymaxcdt1]);
        set(gca,'YTick',ymincdt1:ystepcdt1:ymaxcdt1);
        set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
        if markable
            plot(cdt1timeactual(i,markable),cdt1smooth(cdt1timeactual(i,markable)),'ko','markerfacecolor','k');
        end
    end
end
%%% Remove bad data and sort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodfirstmitosis = mitosesafterspike(:,1)-drugspike;
goodfirstmitosis = goodfirstmitosis(goodfirstmitosis>0);
goodangieg1 = angieg1(angieg1>0);
goodangiemean = angieprevmean(logical(leavein));
goodg1length = g1length(logical(countg1));
%goodg1length = g1length.*countg1;
goodg1length = goodg1length/framesperhr;
goodsg2length = sg2length(logical(leavein));
goodsg2length = goodsg2length/framesperhr;
%%% View Summary Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summarygraph = 0;
if summarygraph
    plot(outputdelay,goodsg2length,'.','markersize',20); %turn on to visualize
    xlabel('Time from Drug Spike --> Mitosis (Hrs)','Fontsize',20);
    ylabel('Time from mitosis --> Cdt1 destruction (Hrs)','Fontsize',20);
    saveas(gcf,'h:\Downloads\Fig.jpg');
end
cd ..\Analysis; %return to this directory
end

function [boundarycheck] = truemax(cdt1,actualtime,numframes,leftedge)
lowertimebound = actualtime-20;
uppertimebound = actualtime+20;
%%% check left side %%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataedge = lowertimebound<=0;
if dataedge
    lowercheck = 1;
else
    stillrising = actualtime==leftedge && cdt1(max([leftedge-3,1]))>cdt1(leftedge);
    slowfall = cdt1(actualtime)-cdt1(lowertimebound)<0.3;
    traceedge = cdt1(actualtime-5)==-10000;
    lowercheck = stillrising || slowfall || traceedge;
end
%%% check right side %%%%%%%%%%%%%%%%%%%%%%%%%%%
dataedge = uppertimebound>numframes;
if dataedge
    uppercheck = 1;
else
    slowfall = cdt1(actualtime)-cdt1(uppertimebound)<0.3;
    traceedge = cdt1(actualtime+5)==-10000;
    uppercheck = dataedge || slowfall || traceedge;
end
boundarycheck = lowercheck || uppercheck;
end
