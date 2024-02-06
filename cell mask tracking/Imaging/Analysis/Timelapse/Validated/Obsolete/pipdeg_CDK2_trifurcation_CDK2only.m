function Align_Mitosis_CDK2(conditions,datadir)
signaloption=1; %1:CDK2 2:pipdeg
minlengthdaughter=60; %default 30
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
%[degstarts,degends]=getdegstartandend_mother(alltracespipdeg,allmotherstats,minlengthmother,relative);
%G2lengths=allmotherstats(:,2)-degends;
%Slengths=degends-degstarts;
%%% measure minimum CDK2 activity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltracesCdk2,1);
earlyval=ones(numtraces,1)*NaN;
lateval=ones(numtraces,1)*NaN;
for i=1:numtraces
    earlyval(i)=alltracesCdk2(i,alldaughterstats(i,1)+10-1); %Sabrina uses 2hr
    lateval(i)=alltracesCdk2(i,alldaughterstats(i,1)+60-1);
end
%figure,hist(earlyval,100);set(gcf,'color','w','PaperPosition',[0 0 4 3]);saveas(gcf,'h:\Downloads\Fig.jpg');
earlycutoff=0.50; %val 6hrs after mitosis
latecutoff=1;

%Cdk2inc=earlyval>earlycutoff & lateval>latecutoff;
Cdk2inc=earlyval>earlycutoff;

%Cdk2low=earlyval<mincutoff;
%figure,hist(lateval(Cdk2low),30);

%Cdk2lowtransient=earlyval<earlycutoff & lateval<earlycutoff;
Cdk2lowtransient=earlyval<earlycutoff;

%Cdk2lowprolonged=lateval<latecutoff;
%notprolonged=lateval>latecutoff;
%figure,hist(earlyval(notprolonged),100);

%{
%%% compare mother G2 lengths by Cdk2inc vs low %%%%%%%%%%%%%%%%%%%%%%%%%%%
G2inc=G2lengths(Cdk2inc)/framesperhr; fprintf('G2inc org n=%0.0f\n',numel(G2inc));
G2low=G2lengths(Cdk2lowtransient)/framesperhr; fprintf('G2low org n=%0.0f\n',numel(G2low));
numnanG2inc=sum(allmotherstats(Cdk2inc,5)>=maxlengthG2 & isnan(G2inc));
numnanG2low=sum(allmotherstats(Cdk2lowtransient,5)>=maxlengthG2 & isnan(G2low));
fprintf('%0.0f traces too long for G2inc measurement\n',numnanG2inc);
fprintf('%0.0f traces too long for G2low measurement\n',numnanG2low);
G2inc(allmotherstats(Cdk2inc,5)<maxlengthG2 | isnan(G2inc))=[];
G2low(allmotherstats(Cdk2lowtransient,5)<maxlengthG2 | isnan(G2low))=[];
G2pval=ranksum(G2inc,G2low);
boxplotdata=[G2inc;G2low];
G2numinc=length(G2inc); G2numlow=length(G2low);
meanG2inc=round(mean(G2inc)*10)/10;
meanG2low=round(mean(G2low)*10)/10;
stringinc=['CDK2inc: mean=',num2str(meanG2inc),'hrs  (n=',num2str(G2numinc),')'];
stringlow=['CDK2low: mean=',num2str(meanG2low),'hrs  (n=',num2str(G2numlow),')'];
offset=[ones(G2numinc,1);2*ones(G2numlow,1)];
figure,boxplot(axes,boxplotdata,offset,'labels',{stringinc,stringlow});
title(['G2 length comparison: p-val=',num2str(G2pval)]);
set(gcf,'color','w','PaperPosition',[0 0 6 4]);
saveas(gcf,'h:\Downloads\Fig_G2lengths.jpg');
%%% compare mother S-phase lengths by Cdk2inc vs low %%%%%%%%%%%%%%%%%%%%%%
Sinc=Slengths(Cdk2inc)/framesperhr; fprintf('Sinc org n=%0.0f\n',numel(Sinc));
Slow=Slengths(Cdk2lowtransient)/framesperhr; fprintf('Slow org n=%0.0f\n',numel(Slow));
numnanSinc=sum(allmotherstats(Cdk2inc,5)>=maxlengthSG2 & isnan(Sinc));
numnanSlow=sum(allmotherstats(Cdk2lowtransient,5)>=maxlengthSG2 & isnan(Slow));
fprintf('%0.0f traces too long for Sinc measurement\n',numnanSinc);
fprintf('%0.0f traces too long for Slow measurement\n',numnanSlow);
Sinc(allmotherstats(Cdk2inc,5)<maxlengthSG2 | isnan(Sinc))=[];
Slow(allmotherstats(Cdk2lowtransient,5)<maxlengthSG2 | isnan(Slow))=[];
Spval=ranksum(Sinc,Slow);
boxplotdata=[Sinc;Slow];
Snuminc=length(Sinc); Snumlow=length(Slow);
meanSinc=round(mean(Sinc)*10)/10;
meanSlow=round(mean(Slow)*10)/10;
stringinc=['CDK2inc: mean=',num2str(meanSinc),'hrs  (n=',num2str(Snuminc),')'];
stringlow=['CDK2low: mean=',num2str(meanSlow),'hrs  (n=',num2str(Snumlow),')'];
offset=[ones(Snuminc,1);2*ones(Snumlow,1)];
figure,boxplot(axes,boxplotdata,offset,'labels',{stringinc,stringlow});
title(['S length comparison: p-val=',num2str(Spval)]);
set(gcf,'color','w','PaperPosition',[0 0 6 4]);
saveas(gcf,'h:\Downloads\Fig_Slengths.jpg');
%}
%%% align daughter data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(alltracesCdk2,1);
dataCdk2daughter=ones(numtraces,numframes)*NaN;
datapipdegdaughter=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    dataCdk2daughter(i,1:alldaughterstats(i,3))=alltracesCdk2(i,alldaughterstats(i,1):alldaughterstats(i,2));
    datapipdegdaughter(i,1:alldaughterstats(i,3))=alltracespipdeg(i,alldaughterstats(i,1):alldaughterstats(i,2));
end
maxlength=prctile(alldaughterstats(:,3),90);
dataCdk2daughter=dataCdk2daughter(:,1:maxlength);
datapipdegdaughter=datapipdegdaughter(:,1:maxlength);
%%% align mother data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataCdk2mother=ones(numtraces,numframes)*NaN;
datapipdegmother=ones(numtraces,numframes)*NaN;
for i=1:numtraces
    dataCdk2mother(i,numframes-allmotherstats(i,3)+1:numframes)=alltracesCdk2(i,allmotherstats(i,1):allmotherstats(i,2));
    datapipdegmother(i,numframes-allmotherstats(i,3)+1:numframes)=alltracespipdeg(i,allmotherstats(i,1):allmotherstats(i,2));
end
mothermax=80;
dataCdk2mother=dataCdk2mother(:,numframes-mothermax+1:numframes);
datapipdegmother=datapipdegmother(:,numframes-mothermax+1:numframes);
%%% concatenate mother & daughter data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if signaloption==1
    datatotal=[dataCdk2mother dataCdk2daughter];
elseif signaloption==2
    datatotal=[datapipdegmother datapipdegdaughter];
end
time=(-mothermax+1:maxlength)/framesperhr;
%%% generate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataCdk2inc=datatotal(Cdk2inc,:);
dataCdk2low=datatotal(Cdk2lowtransient,:);
traceoption=1;
if traceoption
    figure, hold on;
    for i=1:size(dataCdk2inc,1)
        line(time,dataCdk2inc(i,:),'color','b');
    end
    for i=1:size(dataCdk2low,1)
        line(time,dataCdk2low(i,:),'color','r');
    end
    ylim([0 2.5]);
    xlabel('Time relative to Mitosis (hrs)');
    ylabel('CDK2 activity'); ylabel('pipdeg activity');
    set(gcf,'color','w','PaperPosition',[0 0 8 6]);
    saveas(gcf,'h:\Downloads\Fig1.jpg');
    
    figure, hold on;
    binstep=1;
    bincurveshade(time,dataCdk2inc,binstep,'b');
    bincurveshade(time,dataCdk2low,binstep,'r');
    ylim([0 2.5]);
    xlabel('Time relative to Mitosis (hrs)');
    ylabel('CDK2 activity'); ylabel('pipdeg activity');
    set(gcf,'color','w','PaperPosition',[0 0 4 3]);
    saveas(gcf,'h:\Downloads\Fig2.jpg');
end
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
%%% append mother data and record track mother stats %%%%%%%%%%%%%%%%%%%%%%
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
%%% smooth data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numtraces=size(tracedata,1);
DHBnuc=tracedata(:,:,5);
tracesCdk2=tracedata(:,:,7)./tracedata(:,:,5);
tracespipdeg=tracedata(:,:,6);
for i=1:numtraces
    realframes=find(~isnan(DHBnuc(i,:)));
    DHBnuc(i,realframes)=smooth(DHBnuc(i,realframes));
    tracesCdk2(i,realframes)=smooth(tracesCdk2(i,realframes));
    tracespipdeg(i,realframes)=smooth(tracespipdeg(i,realframes));
end
%%% gate traces based on DHB sensor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%badtracespipdeg=gate_pipdeg(tracespipdeg,motherstats,minlengthmother,maxthresh,minthresh);
%badtraces=badtracesCdk2 | badtracespipdeg;
badtraces=badtracesCdk2;
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