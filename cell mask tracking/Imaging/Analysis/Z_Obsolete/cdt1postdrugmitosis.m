%function [outputdata,noentrycount,outputdelay] = cdt1postdrugmitosis(row,col)
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
nocall = 1;
if nocall
    row = [3];
    col = [2];
end
site = [0 1];
framesperhr = 5;
drugspike = 0;          %set to zero if no drug
postdrug1predrug0 = 1;   %this code is written for postdrug mitoses (1)
longestdelay = 5*6; %timescape experiments: 5x6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[signal1,signal2,mitoses] = combinedata(row,col,site,path,drugspike);
tracelist = 1:size(signal2,1);
%%%%%%% remove cells that never mitose %%%%%%%%%%%%%%%%
%fprintf(['count after removing duplicates: ',num2str(length(tracelist)),'\n']);
tracelist(ismember(tracelist,find(max(mitoses,[],2)==0))) = [];
samplesize = length(tracelist);
cdt1 = signal2(tracelist,:);
numframes = size(cdt1,2);
drugtolastframe = numframes-drugspike;
longestentry = drugtolastframe-longestdelay;
cdt1time = zeros(samplesize,1);
actualtime = zeros(samplesize,1);
goodmitoses = mitoses(tracelist,:);
%fprintf(['count after removing non-mitosing cells: ',num2str(length(tracelist)),'\n']);

%%%% mark mitoses and find first mitosis after drug spike %%%%%
mitosisafterspike = zeros(samplesize,1);
nextmitosis = zeros(samplesize,1);
leaveout = zeros(samplesize,1);
noentry = zeros(samplesize,1);
for i=1:samplesize
    mitosistimes = goodmitoses(i,:);
    mitosistimes = mitosistimes(find(mitosistimes>0));
    mitosistimes = sort(mitosistimes);                  %get time window between drugspike & next mitosis
    if postdrug1predrug0
        temp = find(mitosistimes>=drugspike,1);
    else
        temp = find(mitosistimes<drugspike,1,'last');   %only code change required to align by pre-spike mitosis
    end
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

leaveout(mitosisafterspike==0) = 1;
%fprintf(['count after removing non-post-mitosing cells: ',num2str(length(leaveout)-sum(leaveout)),'\n']);
for i=1:samplesize
    timemax = numframes;
    if nextmitosis(i)
        timemax = nextmitosis(i)-3;
    end
    [maxdummy,cdt1time(i)] = max(cdt1(i,mitosisafterspike(i)+1:timemax));
    actualtime(i) = mitosisafterspike(i)+cdt1time(i);
    lowercheck = 0; uppercheck = 0;
    lowertimebound = actualtime(i)-20;
    uppertimebound = actualtime(i)+20;
    if lowertimebound > 0
        if cdt1(i,actualtime(i))-cdt1(i,lowertimebound) < 0.3
            lowercheck = 1;
        end
    else
        lowercheck = 1;
    end
    if uppertimebound <= numframes
        if cdt1(i,actualtime(i))-cdt1(i,uppertimebound) < 0.3
            uppercheck = 1;
        end
    else
        uppercheck = 1;
    end
    boundarycheck = lowercheck * uppercheck;
    if  actualtime(i)+5<=numframes && cdt1(i,actualtime(i)+5)<0
        leaveout(i)=1;
    elseif boundarycheck || actualtime(i)>=numframes-5
        noentry(i)=1;
    end
end
%%%% gate out cells that fall outside our time limits (to avoid biasing) %%%%%
leaveout(mitosisafterspike>(drugspike+longestdelay)) = 1;
leaveout(cdt1time>longestentry)=1;
goodindex = find(leaveout==0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf(['count after removing excess times & neg values: ',num2str(sum(goodindex)),'\n']);
mitosisafterspike = mitosisafterspike-drugspike+1;
cdt1time = (cdt1time)/framesperhr;

outputdata = cdt1time(goodindex);
outputdelay = mitosisafterspike(goodindex)/framesperhr;   %convert to hours
[outputdelay,ix] = sort(outputdelay);     %sort in ascending order
outputdata = outputdata(ix);          %sort delay accordingly
noentrycount = sum(noentry(leaveout==0));
visualize = 0;
if visualize
    plot(outputdelay,outputdata,'.','markersize',20); %turn on to visualize
    xlabel('Time from Drug Spike --> Mitosis (Hrs)','Fontsize',20);
    ylabel('Time from mitosis --> Cdt1 destruction (Hrs)','Fontsize',20);
    %ylim([0.4 1.6]);
    saveas(gcf,'h:\Downloads\Fig.jpg');
end
panelvisualize = 0;
if panelvisualize
    for i=1:length(goodindex)
        cellid = goodindex(i);
        figure(ceil(i/24));             %only 24 plots per figure
        subplot(4,6,mod(i-1,24)+1);
        plot((1:size(cdt1,2)),cdt1(cellid,:),'b');
        if drugspike
            line([drugspike drugspike],[0 2],'Color','r');  %add vertical line when drug is added
        end
        xlim([1 210]); ylim([0 2]);
        hold on;
        plot(actualtime(cellid),cdt1(cellid,actualtime(cellid)),'ro','markerfacecolor','r');
        plot((1:size(angie,2)),angie(cellid,:),'--g');
        plot(mitosisafterspike(cellid)+drugspike,angie(cellid,mitosisafterspike(cellid)+drugspike),'ko','markerfacecolor','k');
    end
end
cd ..\Analysis; %return to this directory