function [outputdata,outputdelay] = postdrugmitosis(row,col,timeafter)
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
%%%%%%%%%%%% Experiment Settings %%%%%%%%
nocall = 0;
if nocall
    row = [3];
    col = [9];
    timeafter = [21 40];
end
site = [0 1];
numframes = 208;
framesperhr = 5;
drugspike = 90;          %set to zero if no drug
postdrug1predrug0 = 1;   %this code is written for postdrug mitoses (1)
%%%%%%%%%%%% Other Settings %%%%%%%%%%%%%
markernum = 3;
tracecolor = 'b';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[traces,mitoses] = combinedata(row,col,site,path,drugspike);
tracelist = 1:size(traces,1);
%%%%%%% remove cells that never mitose %%%%%%%%%%%%%%%%
%fprintf(['count after removing duplicates: ',num2str(length(tracelist)),'\n']);
tracelist(ismember(tracelist,find(max(mitoses,[],2)==0))) = [];
gooddata = traces(tracelist,:);
goodvalue = zeros(size(gooddata,1),1);
goodmitoses = mitoses(tracelist,:);
%fprintf(['count after removing non-mitosing cells: ',num2str(length(tracelist)),'\n']);

%%%% mark mitoses and find first mitosis after drug spike %%%%%
mitosisafterspike = zeros(length(tracelist),1);
nextmitosis = zeros(length(tracelist),1);
leaveout = zeros(length(tracelist),1);
for i=1:size(goodmitoses,1)
    mitosistimes = goodmitoses(i,:);
    for j=1:length(mitosistimes)                    %correct for drugspike offset
        if mitosistimes(j)>=90
            mitosistimes(j) = mitosistimes(j)+1;
        end
    end
    mitosistimes = mitosistimes(find(mitosistimes>0));
    mitosistimes = sort(mitosistimes);                  %get time window between drugspike & next mitosis
    if postdrug1predrug0
        temp = find(mitosistimes>=drugspike+1,1);
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

timepointmin = mitosisafterspike + timeafter(1);
timepointmax = mitosisafterspike + timeafter(end);
leaveout(mitosisafterspike==0) = 1;
%fprintf(['count after removing non-post-mitosing cells: ',num2str(length(leaveout)-sum(leaveout)),'\n']);
for i=1:size(gooddata,1)
    if timepointmax(i)<=208
        if nextmitosis(i) & timepointmax(i)>=nextmitosis(i)
            timepointmax(i)=nextmitosis(i)-2;
        end
        goodvalue(i)=mean(gooddata(i,timepointmin(i):timepointmax(i)));
        if goodvalue(i)<0
            leaveout(i)=1;
        end
    else
        leaveout(i)=1;
    end
end
goodindex = leaveout==0;
%fprintf(['count after removing excess times & neg values: ',num2str(sum(goodindex)),'\n']);
mitosisafterspike = (mitosisafterspike-90)/5;
outputdata = goodvalue(goodindex);
outputdelay = mitosisafterspike(goodindex);
[outputdelay,ix] = sort(outputdelay);     %sort in ascending order
outputdata = outputdata(ix);          %sort delay accordingly
visualize = 0;
if visualize
    plot(outputdelay,outputdata,'.','markersize',20); %turn on to visualize
    xlabel('Time from Drug Spike --> Mitosis (Hrs)','Fontsize',20);
    ylabel('CyclinA-CDK2 Activity','Fontsize',20);
    ylim([0.4 1.6]);
    saveas(gcf,'h:\Downloads\Fig.jpg');
end
cd ..\Analysis; %return to this directory