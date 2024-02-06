clear; close all;
path = 'h:\Documents\Timescape\20120807_12drugs\';
datadir = [path,'Data\'];
%%%%%%%%%%%% Experiment Settings %%%%%%%%
moviename = '3_9_1';
numframes = 208;
framesperhr = 5;
drugspike = 90;          %set to zero if no drug
%%%%%%%%%%%% Other Settings %%%%%%%%%
markernum = 3;

load([datadir,moviename,'_alldata'],'bestsp','best_rc');

%%%%%% store signal values in angieratio and cdt1 matrices
totalcells = size(bestsp{end},1);
totalframes = size(bestsp,3);
angieratio = -10000*ones(totalcells,totalframes);
cdt1 = -10000*ones(totalcells,totalframes);
for f = 1:totalframes
    tempcell = find(bestsp{f}(:,1)~=0);                                     %ignore any cells that have negative coordinates (I'm not sure how negative coordinates would ever appear anyways...)
    angieratio(tempcell,f) = bestsp{f}(tempcell,7)./bestsp{f}(tempcell,5);  %cytosol/nuclear ratio for DHB (Angie's) sensor
    cdt1(tempcell,f) = bestsp{f}(tempcell,6);                               %median nuclear Cdt1 value
end

%%%%%%%%% complete traces %%%%%%%%%
load([path,'GoodCells\',moviename,'_goodcells'],'tracelist');
goodmitoses = zeros(totalcells,10);
generation = zeros(totalcells,1);
for cc = 1:length(tracelist)
    i = tracelist(cc);
    %%% inherit all ancestor cell history and note each mitosis
    temptrace = i;
    m = 1;
    while best_rc(temptrace,5)~=best_rc(temptrace,2)
        mother = best_rc(temptrace,2);
        daughterfirst = best_rc(temptrace,1);
        motherfirst = best_rc(mother,1);
        goodmitoses(i,m) = daughterfirst;
        m = m+1;
        if daughterfirst > drugspike
            generation(i) = generation(i)+1;
        end
        angieratio(i,daughterfirst) = mean([angieratio(i,daughterfirst) angieratio(mother,daughterfirst)]);
        angieratio(i,motherfirst:daughterfirst-1) = angieratio(mother,motherfirst:daughterfirst-1);
        temptrace = mother;
    end
    %%% find each direct descendent mitosis
    daughtermitoses = best_rc(find(best_rc(:,2)==i),1);
    goodmitoses(i,m:m+length(daughtermitoses)-1) = daughtermitoses;
end

%%%%%%%% gate out cells %%%%%%%%%%%%%%%%%%%%
% remove cells that are born two mitoses after drug spike
tracelist(ismember(tracelist,find(generation>1))) = [];
% remove cells that never mitose
tracelist(ismember(tracelist,find(max(goodmitoses,[],2)==0))) = [];
gooddata = angieratio(tracelist,:);
goodmitoses = goodmitoses(tracelist,:);

%%%%% heatmap of angie vs time %%%%%%
drugincluded = [gooddata(:,1:drugspike-1),zeros(length(tracelist),1),gooddata(:,drugspike:totalframes)];
mindata = min(min(drugincluded(find(drugincluded>0))));
offsetdata = drugincluded - mindata;          %for better color display in HeatMap
maxdata = quantile(max(offsetdata,[],2),.85);
specmarkeroffset = maxdata*(markernum/(64-markernum));
offsetdata = offsetdata + specmarkeroffset;
maxdata = maxdata + specmarkeroffset;

cmap = colormap('gray');
cmap(1,:) = [0 1 0];                                %changes last value to green
cmap(2,:) = [1 0 0];                                %changes second to last value to red
cmap(3,:) = [1 1 0];                                %changes third to last value to yellow

drugvalue = maxdata*(1.5)/64;
mitosisvalue = maxdata*(2.5)/64;

offsetdata(:,drugspike) = drugvalue;

spike2mitosis = zeros(length(tracelist),1);
nopostsplit = zeros(length(tracelist),1);
for i=1:size(goodmitoses,1)
    mitosistimes = goodmitoses(i,:);
    for j=1:length(mitosistimes)                    %correct for drugspike offset
        if mitosistimes(j)>=90
            mitosistimes(j) = mitosistimes(j)+1;
        end
    end
    mitosistimes = mitosistimes(find(mitosistimes>0));
    offsetdata(i,mitosistimes) = mitosisvalue;          %mark mitoses green
    mitosistimes = sort(mitosistimes);                  %get time window between drugspike & next mitosis
    temp = find(mitosistimes>=drugspike+1,1);
    if temp
        spike2mitosis(i) = mitosistimes(temp);          %record first mitosis after drugspike
    else
        spike2mitosis(i) = 0;
    end
end

%%%% sort by first mitosis after drug spike %%%%%
[spike2mitosis,ix] = sort(spike2mitosis,'descend');
ordereddata = zeros(length(spike2mitosis),totalframes+1);   %drugspike inserted in as one extra frame
for i=1:length(ix)
    ordereddata(i,:) = offsetdata(ix(i),:);
end

%%%% gate out cells that don't split after drug spike
lastcell = find(spike2mitosis==0,1) - 1;
ordereddata = ordereddata(1:lastcell,:);
spike2mitosis = spike2mitosis(1:lastcell);

totalshift = spike2mitosis(1)-spike2mitosis(end);
aligneddata = zeros(lastcell,size(ordereddata,2)+totalshift);
%%%% align by first mitosis after drug spike %%%%%
for i=1:size(aligneddata,1)
    frontshift = spike2mitosis(1)-spike2mitosis(i);
    endshift = spike2mitosis(i)-(drugspike+1);
    aligneddata(i,:) = [-10000*ones(1,frontshift),ordereddata(i,:),-10000*ones(1,endshift)];
end 

%%%% visualize data%%%%%%%%%%%%%%%
imagesc(aligneddata,[0,maxdata]);
colormap(cmap)


