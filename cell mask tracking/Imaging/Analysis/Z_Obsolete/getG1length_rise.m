function [mitosistime,g1length,g1riselength,quiescent] = getG1length_rise(row,col,g1max,drugspike)
%row=1;col=11;g1max=100;drugspike=40;
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
%    - relevant signal (e.g. angie sensor, Cdt1, etc.)
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
y2min = 0; y2max = 1; y2step = 0.1;
signalwidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[signal1,signal2,mitoses,firstunique] = combinedata_realbirths_history(row,col,site,path);
tracelist = 1:size(signal1,1);  %was signal2
tracelist(ismember(tracelist,find(max(mitoses,[],2)==0))) = [];  %remove any traces with no mitoses
angie = signal1(tracelist,:);
cdt1 = signal2(tracelist,:);
firstunique = firstunique(tracelist,:);
samplesize = length(tracelist);
numframes = size(signal1,2);
xtime = 1:numframes;
maxmitosis = numframes-g1max;  %latest possible mitosis
goodmitoses = mitoses(tracelist,:);
bufsize = 10; tempzeros = zeros(samplesize,bufsize);
mitosesafterspike = tempzeros;    %1st mitosis after drugspike
boundaryflag = tempzeros;
quiescent = tempzeros;
g1length = tempzeros;
sphasestart = tempzeros;
threshpercent = 0.15;
relrisepoint = tempzeros;
actualrisepoint = tempzeros;
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
        smoothcdt1 = smooth(cdt1(i,:));
        smoothcdt1g1 = smoothcdt1(curpoint+1:nextpoint);
        [~,g1length(i,k)] = max(smoothcdt1g1);
        sphasetime = curpoint+g1length(i,k);
        sphasestart(i,k) = sphasetime;
        smoothangie = smooth(angie(i,:));
        
        lowercheck = sphasetime<=curpoint+2;
        uppercheck = sphasetime>=nextpoint-20;
        if uppercheck==0
            nextmin = min(smoothcdt1(sphasetime:nextpoint));
            if nextmin>smoothcdt1(sphasetime)*0.5 && smoothangie(sphasetime)<0.75
                uppercheck=1;
            end
        end
        boundaryflag(i,k) = lowercheck || uppercheck;
        quiescent(i,k) = boundaryflag(i,k) && nextpoint==numframes;
        boundaryflag(i,k) = boundaryflag(i,k) && nextpoint<numframes;
        
        %boundaryflag(i,k) = truemax_postdrop(sphasestart(i,k),curpoint,nextpoint,smoothcdt1);   %check if likely error
        if boundaryflag(i,k)==0
            %g1stretch = smooth(angie(i,curpoint+1:curpoint+g1length(i,k)));
            %g1stretch = smoothangie(curpoint+1:curpoint+g1length(i,k));
            
            %{
            lowg1thresh = min(g1stretch)+range(g1stretch)*threshpercent;
            lowg1 = g1stretch<=lowg1thresh;

            relpos = find(lowg1,1,'first'):find(lowg1,1,'last');
            if maxprev<length(relpos)
                relpos = relpos(end-maxprev+1:end);
            end
            %}
            
            g1slope6 = getslope(smoothangie',6);
            g1slope12 = getslope(smoothangie',12);
            g1slope21 = getslope(smoothangie',21);
            g1depth = 0.5-smoothangie';
            g1dist = 1:length(smoothangie);
            %g1curve = getcurvature(smoothangie',6);
            %risefilter = smooth(g1slope*4+g1curve*2+g1depth*0.5);
            risefilter = smooth(g1slope6*2+g1slope12+g1slope21+g1depth*1+g1dist*0.003);

            relrisepoint(i,k) = find(risefilter(curpoint+1:sphasetime)==max(risefilter(curpoint+1:sphasetime)));
            actualrisepoint(i,k) = curpoint+relrisepoint(i,k);
        end
        k = k+1;
    end
end

%%%% gate out cells that fall outside our time limits (to avoid biasing) %%%%%%%%%%%%%%%%%%%%%%
leavein(mitosesafterspike==0)=0;
leavein(mitosesafterspike>=maxmitosis)=0;
leavein(logical(boundaryflag))=0;  %ignore possible long g1 if cdt1 error

%%%% gate out quiescent traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quiescent(g1length>g1max)=1;        %ignore risetime measurements
%quiescent(leavein==0)=0;           %unnecessary if only keeping leavein
%leavein(quiescent==1)=0;           %only use if leaving out any quiescent

%%%% gate out noisy hDHB traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noisethresh = 30;
noisiness=mean(abs(diff(angie,1,2)),2)*100;
leavein(noisiness>noisethresh,:)=0;

%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 0;
if panelvisualize
    sample=find(sum(leavein,2));
    %for i=1:samplesize
    %sample=[34;61;65;71;110;262];
    for cc=1:length(sample)
        i=sample(cc);
        figure(ceil(cc/20));             %only 24 plots per figure
        set(gcf,'color','w');
        subaxis(4,5,mod(cc-1,20)+1,'ML',0.02,'MR',0.02,'MT',0.03,'MB',0.03,'SH',0.03); %5x4
        signal1trace = angie(i,:);
        signal2trace = 0;
        %signal2trace = smooth(cdt1(i,:));
        
        %%%% uncomment to visualize filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        smoothsignal = smooth(signal1trace);
        g1slope6 = getslope(smoothsignal',6);
        g1slope12 = getslope(smoothsignal',12);
        g1slope21 = getslope(smoothsignal',21);
        g1curve = getcurvature(smoothsignal',12);
        g1depth = 0.5-smoothsignal';
        g1dist = 1:length(smoothsignal);
        %signal2trace = smooth(g1slope*4+g1curve*2+g1depth*0.5+g1dist*0.001);
        signal2trace = smooth(g1slope6*2+g1slope12+g1slope21+g1depth*1+g1dist*0.003);
        %}
        
        [haxes,hline1,hline2] = plotyy(xtime,signal1trace,xtime,signal2trace);
        axes(haxes(1)); hold on;
        axis([0 numframes y1min y1max]);
        set(gca,'YTick',y1min:y1step:y1max);
        set(hline1,'color','b','linewidth',signalwidth);
        title(['Cell ',num2str(i)]);
        
        if drugspike
            line([drugspike drugspike],[y1min y1max],'Color','r');  %add vertical line when drug is added
        end
        scatter(firstunique(i),signal1trace(firstunique(i)),'ko','markerfacecolor','k');
        markable = find(leavein(i,:));
        if markable
            Anaphases = mitosesafterspike(i,markable);
            plot(Anaphases,signal1trace(Anaphases),'ro','markerfacecolor','r');
            scatter(actualrisepoint(i,markable),signal1trace(actualrisepoint(i,markable)),'go','markerfacecolor','g');
        end

        axes(haxes(2)); hold on;
        axis([0 numframes y2min y2max]);
        set(gca,'YTick',y2min:y2step:y2max);
        set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
        line([xtime(1) xtime(end)],[0 0],'Color','k');
    end
end

%%%%%%%%%%%%%%%%%% Compile all good data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g1riselength = g1length-relrisepoint;        %length of g1 rise
mitosistime = mitosesafterspike(logical(leavein));
g1length = g1length(logical(leavein));       %goodg1length = g1length.*countg1;
g1riselength = g1riselength(logical(leavein));
quiescent = quiescent(logical(leavein));
%quiescentcount = sum(sum(quiescent));

%set(gcf,'color','w','PaperPosition',[0 0 10 7]); %3x4 or mini
%saveas(gcf,'h:\Downloads\Fig.jpg');

cd ..\Analysis; %return to this directory