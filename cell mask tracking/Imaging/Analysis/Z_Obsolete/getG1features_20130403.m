%function [mitosistime,g1length,g1riselength,g1toolong] = getG1features(row,col,measurelimit,drugspike)
row=1;col=11;measurelimit=100;drugspike=40;
%%%% Plot signal after post-drugspike mitosis vs drug-mitosis delay %%%%%%%%%%%%%%%%%
% USAGE: - row and col: specify well
%        - measurelimit: max length (in frames) for g1 or g1 latency
%        - drugspike: time (in frames) of drug addition
%        - featurechoice: 1 = G1 length; 2 = G1 latency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Directory Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ..\Functions; %change directory for function calls
%path = 'h:\Documents\Timescape\20120807_12drugs\';
path = 'h:\Documents\Timelapse\Timescape\Steve&Sabrina\';
%%%%%%%%%%%% Experiment Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site = [1];
%%%%%%%%%%%%%% viewer settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1min = 0; y1max = 2.5; y1step = 0.5;
y2min = 0; y2max = 1; y2step = 0.2;
signalwidth = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[signal1,signal2,mitoses,firstunique] = combinedata_realbirths_history(row,col,site,path);
tracelist = 1:size(signal1,1);  %was signal2
tracelist(ismember(tracelist,find(max(mitoses,[],2)==0))) = [];  %remove any traces with no mitoses
angie = signal1(tracelist,:);
cdt1 = signal2(tracelist,:);
firstunique = firstunique(tracelist,:);
samplesize = length(tracelist);
numframes = size(signal1,2);
xtime = 1:numframes;
measuremargin = 10;   %previously 20
maxmitosis = numframes-measurelimit-measuremargin;  %latest possible mitosis
goodmitoses = mitoses(tracelist,:);
tempzeros = zeros(samplesize,10);
mitosesafterspike = tempzeros;    %1st mitosis after drugspike
Sflag = tempzeros;
G1flag = tempzeros;
g1toolong = tempzeros;
g0toolong = tempzeros;
g1length = tempzeros;
sphaseentry = tempzeros;
threshpercent = 0.15;
g0length = tempzeros;
g1entry = tempzeros;
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
        sphaseentry(i,k) = sphasetime;    %only for visualization
        smoothangie = smooth(angie(i,:));
        lowercheck = sphasetime<=curpoint+2;
        uppercheck = sphasetime>=nextpoint-measuremargin;
        if uppercheck==0
            %nextmin = min(smoothcdt1(sphasetime:nextpoint));
            %if nextmin>smoothcdt1(sphasetime)*0.5 && smoothangie(sphasetime)<0.75
            if smoothangie(sphasetime)<0.75
                uppercheck=1;
            end
        end
        Sflag(i,k) = lowercheck || uppercheck;
        g1toolong(i,k) = Sflag(i,k) && nextpoint==numframes;
        Sflag(i,k) = Sflag(i,k) && nextpoint<numframes;
        
        if Sflag(i,k)==0
            g1slope6 = getslope(smoothangie',6); g1slope12 = getslope(smoothangie',12); g1slope20 = getslope(smoothangie',20);
            %g1slope = g1slope6*2+g1slope12+g1slope20;
            %g1slope = getslope_array(smoothangie',[6,6,12,20]);
            %g1slope = getslope_array(smoothangie',[6:10,6:10,11:20]);
            g1slope = getslope_array(smoothangie',[1:10]);
            %g1slope = getslope_array(smoothangie',[1:20]);
            g1depth = 0.5-smoothangie';
            g1dist = 1:length(smoothangie);
            risefilter = smooth(g1slope*0.5+g1depth+g1dist*0.003);
            if g1toolong(i,k)==1
                risefilter = smooth(g1slope+g1depth+g1dist*0.003);   %more sensitive to gradient
                sphasetime = numframes;
            end
            %risefilter = smooth(g1slope*0.5+g1depth+g1dist*0.003);
            %risefilter = smooth(g1slope*2+g1depth*5+g1dist*0.01);
            g0length(i,k) = find(risefilter(curpoint+1:sphasetime)==max(risefilter(curpoint+1:sphasetime)),1,'first');
            g1entry(i,k) = curpoint+g0length(i,k);   %for visualization only
            
            %G1flag(i,k) = something;
        end
        k = k+1;
    end
end

%%%% gate out cells that fall outside our time limits (to avoid biasing) %%%%%%%%%%%%%%%%%%%%%%
leavein(mitosesafterspike==0)=0;
leavein(mitosesafterspike>=maxmitosis)=0;
leavein(logical(Sflag))=0;  %errant trace
%leavein(g1toolong==0)=0;   %toggle to visualize subset

%%%% gate out noisy hDHB traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noisethresh = 30;
noisiness=mean(abs(diff(angie,1,2)),2)*100;
leavein(noisiness>noisethresh,:)=0;

%%%% label quiescent traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g1toolong(g1length>measurelimit)=1;        %ignore risetime measurements
g0toolong(g0length>measurelimit)=1;

%%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
panelvisualize = 1;
if panelvisualize
    sample=find(sum(leavein,2));
    %for i=1:samplesize
    %sample=[5;8;15;34;35;57;65;71;74;110;115;137;143;160;165;194;197;201;210;225;226;231;238;239;243;245;248;262;264;266;269;273];
    for cc=1:length(sample)
        i=sample(cc);
        figure(ceil(cc/20));             %only 24 plots per figure
        set(gcf,'color','w');
        set(gcf,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]);
        subaxis(4,5,mod(cc-1,20)+1,'ML',0.02,'MR',0.02,'MT',0.03,'MB',0.03,'SH',0.03); %5x4
        signal1trace = angie(i,:);
        signal2trace = 0;
        %signal2trace = smooth(cdt1(i,:));
        
        %%%% uncomment to visualize filters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %{
        smoothsignal = smooth(signal1trace);
        %g1slope6 = getslope(smoothsignal',6); g1slope12 = getslope(smoothsignal',12); g1slope20 = getslope(smoothsignal',20);
        %g1slope = getslope_array(smoothsignal',[6:10,6:10,11:20]);
        g1slope = getslope_array(smoothsignal',[1:10]);
        %g1slope = getslope_array(smoothsignal',[6,6,12,20]);
        g1depth = 0.5-smoothsignal';
        g1dist = 1:length(smoothsignal);
        %signal2trace = smooth(g1slope6*2+g1slope12+g1slope20+g1depth+g1dist*0.003);
        signal2trace = smooth(g1slope*0.5+g1depth+g1dist*0.003);
        %signal2trace = smooth(g1slope+g1depth*2);
        
        %signal2trace = smooth(g1slope*0.5+g1depth+g1dist*0.003);
        %signal2trace = smooth(g1slope*2+g1depth*5+g1dist*0.01);
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [haxes,hline1,hline2] = plotyy(xtime,signal1trace,xtime,signal2trace);
        axes(haxes(1)); hold on;
        axis([0 numframes y1min y1max]);
        set(gca,'YTick',y1min:y1step:y1max);
        set(hline1,'color','b','linewidth',signalwidth);
        title(['Cell ',num2str(i)]);
        
        if drugspike
            line([drugspike drugspike],[y1min y1max],'Color','r');  %add vertical line when drug is added
        end
        
        markable = find(leavein(i,:));
        scatter(mitosesafterspike(i,markable),signal1trace(mitosesafterspike(i,markable)),'ro','markerfacecolor','r');
        scatter(g1entry(i,markable),signal1trace(g1entry(i,markable)),'go','markerfacecolor','g');

        scatter(firstunique(i),signal1trace(firstunique(i)),'ko','markerfacecolor','k');
        axes(haxes(2)); hold on;
        %scatter(sphaseentry(i,markable),signal2trace(sphaseentry(i,markable)),'go','markerfacecolor','g');
        axis([0 numframes y2min y2max]);
        set(gca,'YTick',y2min:y2step:y2max);
        set(hline2,'linestyle','--','color','g','linewidth',signalwidth);
    end
end

%%%%%%%%%%%%%%%%%% Compile all good data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g1riselength = g1length-g0length;        %length of g1 rise
mitosistime = mitosesafterspike(logical(leavein));
g1length = g1length(logical(leavein));       %goodg1length = g1length.*countg1;
g1riselength = g1riselength(logical(leavein));
g1toolong = g1toolong(logical(leavein));
%quiescentcount = sum(sum(quiescent));

%set(gcf,'color','w','PaperPosition',[0 0 10 7]); %3x4 or mini
%saveas(gcf,'h:\Downloads\Fig.jpg');

cd ..\Analysis; %return to this directory