function main()
clear; close all;
%%%%%%%%%%%% Experiment Settings %%%%%%%%
numframes = 208;
framesperhr = 5;
drugspike = 90;          %set to zero if no drug
post1pre0 = 1;
cut = 1;
cycletime = 0;
%%%%%%%%%%%% Other Settings %%%%%%%%%
markernum = 3;

row = [3];
col = [9];
site = [0 1];
[traces,mitoses] = combinedata(row,col,site);
tracelist = 1:size(traces,1);
%%%%%%%% gate out cells %%%%%%%%%%%%%%%%%%%%
% remove cells that never mitose
fprintf(['count after removing duplicates: ',num2str(length(tracelist)),'\n']);
tracelist(ismember(tracelist,find(max(mitoses,[],2)==0))) = [];
gooddata = traces(tracelist,:);
goodmitoses = mitoses(tracelist,:);
fprintf(['count after removing non-mitosing cells: ',num2str(length(tracelist)),'\n']);

%%%% mark drugspike and set colors %%%%%%
drugincluded = [gooddata(:,1:drugspike-1),zeros(length(tracelist),1),gooddata(:,drugspike:numframes)];
mindata = min(min(drugincluded(find(drugincluded>0))));
mindata = 0.5;
offsetdata = drugincluded - mindata;          %for better color display in HeatMap
maxdata = quantile(max(offsetdata,[],2),.50);
maxdata = 0.8;
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
    offsetdata(i,mitosistimes) = mitosisvalue;          %mark mitoses green
    mitosistimes = sort(mitosistimes);                  %get time window between drugspike & next mitosis
    if post1pre0
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
ordereddata = zeros(length(mitosis2spike),numframes+1);   %drugspike inserted in as one extra frame
for i=1:length(ix)
    ordereddata(i,:) = offsetdata(ix(i),:);
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
cutofflength = 60;
cutoff = find(endshift<cutofflength,1);
if cut
    aligneddata = aligneddata(cutoff:end,:);
    fprintf(['after removing below cutoff: ',num2str(size(aligneddata,1)),'\n']);

    %%%% ignore first <totalshift> and last <cutoff> frames %%%%
    aligneddata = aligneddata(:,totalshift+1:end-cutofflength);
end

%%%% visualize data%%%%%%%%%%%%%%%
timebefore = -mitosis2spike(end)/5;
timeafter = (size(aligneddata,2)-mitosis2spike(end))/5;
imagesc([timebefore timeafter],1,aligneddata,[0,maxdata]);

colormap(cmap)
if ~cut
    line([1 extendedframes],[cutoff cutoff],'Color','m');
end
set(0,'Units','normalized');                %sets screensize units by pixels
set(gcf,'Units','normalized');
set(gcf,'Position',[0.3 0.3 0.4 0.3]);
set(gca,'Position',[0.05 0.15 0.9 0.8]);
set(gca,'YTick',[]);
xlabel('Relative Time (Hrs)','Fontsize',12);
saveas(gcf,'h:\Downloads\Fig.jpg');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [traces, mitoses] = combinedata(rowmat,colmat,sitemat)
path = 'h:\Documents\Timescape\20120807_12drugs\';
datadir = [path,'Data\'];
gatedir = [path,'GoodCells\'];
drugspike = 90;
traces = [];
mitoses = [];
duplicates = 0;
for row=rowmat
    for col=colmat
        for site=sitemat
            moviename = [num2str(row),'_',num2str(col),'_',num2str(site)];
            load([datadir,moviename,'_alldata'],'best_rc');
            load([gatedir,moviename,'_goodcells'],'tracelist');
            load([gatedir,moviename,'_singletracedata'],'signal1','goodmitoses');
            for cc=1:length(tracelist)
                i=tracelist(cc);
                eachmitosis = goodmitoses(i,:);
                postspikemitoses = sort(eachmitosis(eachmitosis>drugspike));
                if length(postspikemitoses)>1
                    if best_rc(i,1)==postspikemitoses(2) & ismember(best_rc(i,2),tracelist) %'duplicate' data
                        duplicates = duplicates+1;
                        tracelist(cc)=0;
                    end
                end
            end
            tracelist = tracelist(tracelist>0);
            traces = [traces;signal1(tracelist,:)];
            mitoses = [mitoses;goodmitoses(tracelist,:)];
        end
    end
end
fprintf(['duplicate count: ',num2str(duplicates),'\n']);
end

