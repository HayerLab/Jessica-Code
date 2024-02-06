function Hoechst_EdU_Ab(conditions,datadir)
panel=0; %0:single image
condnum=size(conditions,1);
%colorcode=distributecolormap(jet,condnum);
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    Ab1vals=[];
    Ab2vals=[];
    Ab3vals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [Ab1val,Ab2val,Ab3val]=main(datadir,shot);
                Ab1vals=[Ab1vals;Ab1val];
                Ab2vals=[Ab2vals;Ab2val];
                Ab3vals=[Ab3vals;Ab3val];
            end
        end
    end
    valuesx=Ab1vals;
    valuesy=Ab2vals;
    valuesz=Ab3vals;
    %maxz=3.5; minz=2.9; midz=3.20; %pRbAbChar 20131217 BGsub offset 3.5 2.9 3.2
    %maxz=3.2; minz=2.80; midz=2.96; %serumrelease doublerabbit (median) 3.2 2.6 2.96 (2.91)
    %maxz=2.7; minz=2.3; midz=2.5; %serumrelease 20131126 (10perc10nucr median)
    maxz=6.5; minz=4.5; midz=5.75; %serumrelease 20131126 (bgjustmean mean)
    %hist(valuesz,1.95:0.05:6.95);xlim([2 7]);
    %sampleid=find(valuesz>maxz);
    sampleid=find(valuesz>(midz-0.25) & valuesz<(midz+0.25));
    %sampleid=find(valuesz>(midz-0.15) & valuesz<(midz+0.15));
    %sampleid=find(valuesz<minz);
    valuesx(sampleid)=[]; valuesy(sampleid)=[]; valuesz(sampleid)=[];
    valuesz(valuesz>maxz)=maxz; valuesz(valuesz<minz)=minz;

%     figure; hist(Ab3vals,(minz-0.02):0.02:(maxz+0.02));xlim([minz maxz]);
%     set(gcf,'color','w','PaperPosition',[0 0 6 5]);
%     saveas(gcf,'h:\Downloads\Fig_Hist.jpg');
%     figure; imagesc(1:128); colormap(cmap);
%     set(gcf,'color','w','PaperPosition',[0 0 6 1]);
%     saveas(gcf,'h:\Downloads\Fig_Colorbar.jpg');
    
    if panel==1
        subaxis(3,4,i,'ML',0.02,'MR',0.01,'MT',0.2,'MB',0.1,'SH',0.02); %5x4
    end
    %make3dmap(valuesx,valuesy,valuesz,minz,midz,maxz);
    %makedensitymap(valuesx,valuesy);
    makeprobmap(valuesx,valuesy,valuesz,midz);

    hold on; line([3.5 10],[10 3.5],'linewidth',2,'color','r');  %
    axis([4 9 3 9]); %[3 10 3 10]
    title(char(conditions(i,1)),'FontSize',16);
end

%xstring='p21 log2(meanRFU)'; ystring='CycD1 log2(meanRFU)';
%xlabel(xstring); ylabel(ystring);
if panel==0
    paperwidth=6; paperheight=5;
else
    paperwidth=12; paperheight=10;
end
set(gcf,'color','w','PaperPosition',[0 0 paperwidth paperheight]); %6 5 (small panel: 4 3.3)
saveas(gcf,'h:\Downloads\Fig.jpg');
%{
figure,boxplot(axes,boxplotdata,offset,'labels',namecount);
title(titlestring);ylabel(descstring);
ylim([bmin bmax]);
set(gcf,'color','w','PaperPosition',[0 0 4 4]);
saveas(gcf,'h:\Downloads\Fig_boxplots.jpg');
%}
%{
figure;
errorbar(sigtime,sigmean,sigstd,'-mo','linewidth',2,'markeredgecolor','k','markerfacecolor',[.49 1 .63],'markersize',10);
title(titlestring);xlabel('time since serum release (hr)'); ylabel(descstring);
set(gcf,'color','w','PaperPosition',[0 0 4 4]);
saveas(gcf,'h:\Downloads\Fig_boxplots.jpg');
%}
end


function [Ab1,Ab2,Ab3]=main(datadir,shot)
%1:X 2:Y 3:Area 4:Int(Hoechst) 5:med(pRb) 6:med(p21) 7:med(CycD) 8-10:mean 11-13:max
load([datadir,shot,'.mat'],'IFdata');
%Hoechstval=IFdata(:,3).*IFdata(:,8);
Hoechstval=IFdata(:,4);
%Ab1val=IFdata(:,5)-IFdata(:,14);
Ab1val=IFdata(:,6); %perc10nucr10:6 org5 med5 mean9 max13
Ab1val(Ab1val<1)=1; Ab1val=log2(Ab1val);
%Ab2val=IFdata(:,6)-IFdata(:,15);
Ab2val=IFdata(:,7); %p10n10:7 org6 med6 mean10 max14
Ab2val(Ab2val<1)=1; Ab2val=log2(Ab2val);
Ab3val=IFdata(:,8); %p10n10:8 org7 med7 mean11 max15
Ab3val(Ab3val<1)=1; Ab3val=log2(Ab3val);
%%% Gate by Hoechst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(Hoechstval,0:10000:1000000);
G1cells=find(Hoechstval<275000); %pRbAbChar 20131217 Bgsub G1:<250000 G2:400000<x<600000
%G1cells=find(Hoechstval>400000 & Hoechstval<600000);
%G1cells=find(Hoechstval<225000 & Hoechstval>100000); %G1 default 200000
%G1cells=find(Hoechstval>300000); %G2
%G1cells=find(Hoechstval<600000); %All
Ab1=Ab1val(G1cells);
Ab2=Ab2val(G1cells);
Ab3=Ab3val(G1cells);
end

function make3dmap(valuesx,valuesy,valuesz,minz,midz,maxz)
cmapz=mat2gray(valuesz)*128; cmapz=round(cmapz); cmapz(cmapz<1)=1; cmapz(cmapz>128)=128;
lowcode=[0 0 1]; highcode=[1 1 0]; midcode=mean([lowcode;highcode]);
cmapmid=round(128*(midz-minz)/(maxz-minz));
cmap=makecmap_midpoint(lowcode,midcode,highcode,cmapmid);
colorsz=cmap(cmapz,:);
scatter(valuesx,valuesy,30,colorsz,'fill');
end

function makedensitymap(valuesx,valuesy)
stepsize=0.2;
halfstep=0.5*stepsize;
binx=min(valuesx)+halfstep:stepsize:max(valuesx)-halfstep;
biny=min(valuesy)+halfstep:stepsize:max(valuesy)-halfstep;
numbinx=length(binx); numbiny=length(biny);
vbinx=zeros(numbinx*numbiny,1); vbiny=vbinx; vbinz=vbinx;
total=numel(valuesx);
cc=0;
for i=1:numbinx
    samplex=find(valuesx>=binx(i)-halfstep & valuesx<binx(i)+halfstep);
    for j=1:numbiny
        cc=cc+1;
        vbinx(cc)=binx(i);
        vbiny(cc)=biny(j);
        if isempty(samplex)
            continue;
        end
        samplebinidx=valuesy(samplex)>=biny(j)-halfstep & valuesy(samplex)<biny(j)+halfstep;
        samplebin=samplex(samplebinidx);
        samplecount=numel(samplebin);
        if samplecount<=1
            samplecount=0;
        end
        vbinz(cc)=100*samplecount/total;
    end
end
cmapz=mat2gray(vbinz)*64; cmapz=round(cmapz); cmapz(cmapz<1)=1; cmapz(cmapz>64)=64;
cmap=jet;
cmap(1,:)=[1 1 1]; %0:white
colorsz=cmap(cmapz,:);
scatter(vbinx,vbiny,80,colorsz,'fill','square');
%imagesc(countmap);
%cmap=colormap(jet);
%cmap(1,:)=[1 1 1];
%colormap(cmap);
%colorbar;
end

function makeprobmap(valuesx,valuesy,valuesz,midz)
stepsize=0.2;
halfstep=0.5*stepsize;
binx=min(valuesx)+halfstep:stepsize:max(valuesx)-halfstep;
biny=min(valuesy)+halfstep:stepsize:max(valuesy)-halfstep;
numbinx=length(binx); numbiny=length(biny);
vbinx=zeros(numbinx*numbiny,1); vbiny=vbinx; vbinz=vbinx;
cc=0;
for i=1:numbinx
    samplex=find(valuesx>=binx(i)-halfstep & valuesx<binx(i)+halfstep);
    for j=1:numbiny
        cc=cc+1;
        vbinx(cc)=binx(i);
        vbiny(cc)=biny(j);
        if isempty(samplex)
            continue;
        end
        samplebinidx=valuesy(samplex)>=biny(j)-halfstep & valuesy(samplex)<biny(j)+halfstep;
        samplebin=samplex(samplebinidx);
        samplecount=numel(samplebin);
        samplepos=sum(valuesz(samplebin)>midz);
        if samplecount<=10
            continue;
        end
        vbinz(cc)=100*samplepos/samplecount;
    end
end
cmapz=mat2gray(vbinz)*64; cmapz=round(cmapz); cmapz(cmapz<1)=1; cmapz(cmapz>64)=64;
cmap=jet;
cmap(1,:)=[1 1 1]; %0:white
colorsz=cmap(cmapz,:);
scatter(vbinx,vbiny,80,colorsz,'fill','square');
end