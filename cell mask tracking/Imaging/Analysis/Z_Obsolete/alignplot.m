cd ..\Functions; %change directory for function calls
path = 'h:\Documents\Timescape\20120807_12drugs\';
%%%%%%%%%%%% Experiment Settings %%%%%%%%
numframes = 208;
framesperhr = 5;
drugspike = 90;          %set to zero if no drug
postdrug1predrug0 = 1;
cut = 1;
cycletime = 0;
averagemode = 1;
%%%%%%%%%%%% Other Settings %%%%%%%%%
markernum = 3;

row = [3];
col = [2];
site = [0 1];
hold on;
tracecolor = 'b';
[traces1,traces2,mitoses] = combinedata(row,col,site,path,drugspike);
tracelist = 1:size(traces1,1);
%%%%%%%% gate out cells %%%%%%%%%%%%%%%%%%%%
% remove cells that never mitose
fprintf(['count after removing duplicates: ',num2str(length(tracelist)),'\n']);
tracelist(ismember(tracelist,find(max(mitoses,[],2)==0))) = [];
gooddata = traces1(tracelist,:);
goodmitoses = mitoses(tracelist,:);
fprintf(['count after removing non-mitosing cells: ',num2str(length(tracelist)),'\n']);

%%%% mark mitoses and find first mitosis after drug spike %%%%%
mitosis2spike = zeros(length(tracelist),1);
spike2mitosis = zeros(length(tracelist),1);
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
        mitosis2spike(i) = mitosistimes(temp);          %record mitosis
        if temp < length(mitosistimes)                  %if there is a subsequent mitosis (only relevant for pre-spike mitoses)
            spike2mitosis(i) = mitosistimes(temp+1);
        else
            spike2mitosis(i) = 0;
        end
    else
        mitosis2spike(i) = 0;
    end
end

if cycletime
    bothexist = mitosis2spike & spike2mitosis;
    gaptimes = zeros(sum(bothexist),2);
    gaptimes(:,1) = drugspike - mitosis2spike(bothexist);
    gaptimes(:,2) = spike2mitosis(bothexist)-drugspike;
    gaptimes
    save('h:\Downloads\gaptimesLY29.mat','gaptimes');
end

%%%% sort by first mitosis after drug spike %%%%%
[mitosis2spike,ix] = sort(mitosis2spike,'descend');
ordereddata = zeros(length(mitosis2spike),numframes);   %drugspike inserted in as one extra frame
for i=1:length(ix)
    ordereddata(i,:) = gooddata(ix(i),:);
end

%%%% gate out cells that don't split after drug spike
lastcell = find(mitosis2spike==0,1) - 1;
ordereddata = ordereddata(1:lastcell,:);
mitosis2spike = mitosis2spike(1:lastcell);
fprintf(['after removing post-spike non-mitotic cells: ',num2str(length(mitosis2spike)),'\n']);

totalshift = mitosis2spike(1)-mitosis2spike(end);
extendedframes = size(ordereddata,2)+totalshift;
aligneddata = zeros(lastcell,extendedframes);

%%%% align by first mitosis after drug spike %%%%%
endshift = zeros(size(aligneddata,1),1);
for i=1:size(aligneddata,1)
    frontshift = mitosis2spike(1)-mitosis2spike(i);
    endshift(i) = mitosis2spike(i)-mitosis2spike(end);
    aligneddata(i,:) = [-10000*ones(1,frontshift),ordereddata(i,:),-10000*ones(1,endshift(i))];
end

%%%% only look at traces below cutoff %%%%
cutofflength = 100;
cutoff = find(endshift<cutofflength,1);
if cut
    aligneddata = aligneddata(cutoff:end,:);
    fprintf(['after removing below cutoff: ',num2str(size(aligneddata,1)),'\n']);

    %%%% ignore first <totalshift> and last <cutoff> frames %%%%
    aligneddata = aligneddata(:,totalshift+1:end-cutofflength);
end

%%%% visualize data %%%%%%%%%%%%%%%
timebefore = -mitosis2spike(end)/5;
timeafter = (size(aligneddata,2)-mitosis2spike(end)-1)/5;

if averagemode
    line([timebefore:0.2:timeafter],mean(aligneddata),'color',tracecolor);
else
    for i=1:size(aligneddata,1)
        line([timebefore:0.2:timeafter],smooth(aligneddata(i,:)),'color',tracecolor);
    end
end

if ~cut
    line([1 extendedframes],[cutoff cutoff],'Color','m');
end
ylim([0.3 2]);
set(0,'Units','normalized');                %sets screensize units by pixels
set(gcf,'Units','normalized');
set(gcf,'Position',[0.55 0.4 0.4 0.5]);
set(gca,'Position',[0.05 0.15 0.9 0.8]);
xlabel('Relative Time (Hrs)','Fontsize',12);
saveas(gcf,'h:\Downloads\Fig.jpg');

timekeeper = timebefore:0.2:timeafter;
timerange = find(timekeeper<6.1 & timekeeper>4.9);
intensities = aligneddata(:,timerange);
save('h:\Downloads\intensitiespkcb.mat','intensities');
cd ..\Analysis; %return to this directory
