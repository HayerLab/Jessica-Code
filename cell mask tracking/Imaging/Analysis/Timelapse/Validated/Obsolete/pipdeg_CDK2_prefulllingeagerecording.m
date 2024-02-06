function pipdeg_CDK2(conditions,datadir)
minlengthdaughter=30;
minlengthmother=10;
maxlengthG2=20*5;
maxlengthSG2=30*5;
condnum=size(conditions,1);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    alltracesCdk2=[];
    alltracespipdeg=[];
    alldaughterstats=[];
    allmotherstats=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                %shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [tracesCdk2,tracespipdeg,daughterstats,motherstats]=main(datadir,shot,minlengthdaughter,minlengthmother);
                alltracesCdk2=[alltracesCdk2;tracesCdk2];
                alltracespipdeg=[alltracespipdeg;tracespipdeg];
                alldaughterstats=[alldaughterstats;daughterstats];
                allmotherstats=[allmotherstats;motherstats];
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
framesperhr=5;
[numtraces,numframes]=size(alltracesCdk2);
%sampleid=(1:numtraces)';
sampleid=ones(numtraces,1);
%%% extract pipdeg features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relative=0;
[degstarts,degends]=getdegstartandend_smooth(alltracespipdeg,allmotherstats,minlengthmother,relative);
G2lengths=allmotherstats(:,2)-degends;
Slengths=degends-degstarts;
%%% categorize daughter CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltracesCdk2,1);
earlyval=ones(numtraces,1)*NaN;
maxval=ones(numtraces,1)*NaN;
lastval=ones(numtraces,1)*NaN;
for i=1:numtraces
    earlyval(i)=alltracesCdk2(i,alldaughterstats(i,1)+minlengthdaughter-1); %Sabrina uses 2hr
    maxval(i)=max(alltracesCdk2(i,alldaughterstats(i,1)+minlengthdaughter-1:alldaughterstats(i,2))); %max from minlength-end
    lastval(i)=alltracesCdk2(i,alldaughterstats(i,2));
end
%figure,hist(minval,100);
mincutoff=0.60; %minval at 20frames
maxcutoff=0.60; %maxval cutoff F20-end
Cdk2inc=earlyval>mincutoff & lastval>maxcutoff+0.4; %default 0.4
Cdk2low=earlyval<mincutoff & maxval<maxcutoff;
%%% categorize mother phase lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure,hist(G2lengths,100);
G2short=G2lengths>10 & G2lengths<30;
G2long=G2lengths>50 & G2lengths<80;
%%% compare mother phase lengths by daughter CDK2 fate %%%%%%%%%%%%%%%%%%%%
motherG2ofdaughterCDK2fate(G2lengths,Cdk2inc,Cdk2low,allmotherstats,maxlengthG2,framesperhr);
motherSofdaughterCDK2fate(Slengths,Cdk2inc,Cdk2low,allmotherstats,maxlengthSG2,framesperhr);
%%% align traces to mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daughtermax=prctile(alldaughterstats(:,3),90);
mothermax=prctile(allmotherstats(:,3),90);
[time,datatotal]=alignmotherdaughter(alltracesCdk2,alldaughterstats,allmotherstats,daughtermax,mothermax,framesperhr);
%%% plot traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xstring='Time relative to Mitosis (hrs)';
ystring='CDK2 activity'; ylimits=[0 2.5];
comparetrends(datatotal,time,Cdk2inc,Cdk2low,xstring,ystring,ylimits);
%comparetrends(datatotal,time,G2short,G2long,xstring,ystring,ylimits);
%%% report number of traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%0.0f traces\n',numtraces);
end


function [tracesCdk2,tracespipdeg,daughterstats,motherstats]=main(datadir,shot,minlengthdaughter,minlengthmother)
%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','genealogy','jitters');
tracestats=getstats(tracedata,genealogy); %start,end,length,genealogy
%%% get cells that are daughters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplecells=find(~isnan(genealogy));
%%% get full durations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fulltracedata=linkancestry(tracedata,tracestats,samplecells);
%%% append mother data and record mother stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%
motherstats=ones(size(tracestats,1),size(tracestats,2)+1)*NaN;
for i=1:numel(samplecells)
    s=samplecells(i);
    m=tracestats(s,4);
    ancestrylength=tracestats(m,2)-find(~isnan(fulltracedata(s,:,1)),1,'first')+1;
    motherstats(s,:)=[tracestats(m,:) ancestrylength];
    motherframes=motherstats(s,1):motherstats(s,2);
    tracedata(s,motherframes,:)=tracedata(m,motherframes,:);
end
tracedata=tracedata(samplecells,:,:);
daughterstats=tracestats(samplecells,:);
motherstats=motherstats(samplecells,:);
%%% smooth CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
DHBnuc=tracedata(:,:,5);
tracesCdk2=tracedata(:,:,7)./tracedata(:,:,5);
tracespipdeg=tracedata(:,:,6);
for i=1:numtraces
    realframes=find(~isnan(DHBnuc(i,:)));
    DHBnuc(i,realframes)=smooth(DHBnuc(i,realframes));
    tracesCdk2(i,realframes)=smooth(tracesCdk2(i,realframes));
end
%%% gate CDK2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(DHBnuc(:),0:10:2010); xlim([0 2000]);
%hist(max(DHBnuc,[],2),0:20:1020); xlim([0 1000]);
thresh=75; %HW:150 pipdeg0715/0719:localbg=75 sigbg=50
%badtracesCdk2=gate_Cdk2(DHBnuc,tracesCdk2,daughterstats,minlengthdaughter,thresh);
badtracesCdk2=gate_Cdk2_noisecalcdiff1(DHBnuc,tracesCdk2,daughterstats,minlengthdaughter,thresh);
%%% gate traces based on pipdeg sensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(max(tracespipdeg,[],2),0:20:2020); xlim([0 2000]);
%hist(min(tracespipdeg,[],2),0:2:202); xlim([0 200]);
%hist(log2(min(tracespipdeg,[],2)),0:0.1:12.1); xlim([0 12]);
maxthresh=150; %50prc:50 10nuc50prc:100-->200
minthresh=75; %50prc:25 10nuc50prc:100-->75
badtracespipdeg=gate_pipdeg(tracespipdeg,motherstats,minlengthmother,maxthresh,minthresh);
badtraces=badtracesCdk2 | badtracespipdeg;
%%% remove non-sample cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tracesCdk2=tracesCdk2(~badtraces,:);
tracespipdeg=tracespipdeg(~badtraces,:);
daughterstats=daughterstats(~badtraces,:);
motherstats=motherstats(~badtraces,:);
%%% normalize pipdeg traces by max in mother %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(tracespipdeg,1)
    maxval=max(tracespipdeg(i,motherstats(i,1):motherstats(i,2)));
    tracespipdeg(i,:)=tracespipdeg(i,:)/maxval;
end
end