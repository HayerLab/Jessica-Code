function Hoechst_EdU_Ab(conditions,datadir)
condnum=size(conditions,1);
boxplotdata=[];
offset=[];
%%% Mouse log2
%bmin=0;  %G1: med 7  max 7  int 16     S: med:7 max 7 sum 16    G2: med 7  max 7 sum 16
%bmax=11; %G1: med 10 max 12 int 19     S: med:11 max 14 sum 20   G2: med 11 max 14 sum 21
%%% Goat raw
% bmin=50;    %G1: med 50     S: med 0    G2: med 0 max 100
% bmax=450;  %G1: med 350    S: med 400  G2: med 400 max 800
%%% p21
bmin=5; %6 CycD 4
bmax=12; %10 CycD 11
bstep=(bmax-bmin)/50; %pdf 30
bin=bmin:bstep:bmax;
namecount=cell(condnum,1);
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
    valuesx=Ab2vals;
    valuesy=Ab1vals;
    valuesz=Ab3vals;
    %maxz=3.55; minz=2.80; midz=3.20; %pRbAbChar 20131217
    maxz=3.50; minz=3.20; midz=3.35; %pRbAbChar 20131217 BGsub
    %maxz=3.2; minz=2.7; midz=2.97; %serum release 20131216 bgtemp mean
    %maxz=3.35; minz=2.9; midz=3.13; %serum release 20131216 centroid max
    %maxz=12; minz=6; midz=10;
    sampleid=find(valuesz>maxz);
    valuesx(sampleid)=[]; valuesy(sampleid)=[]; valuesz(sampleid)=[];
    valuesz(valuesz>maxz)=maxz; valuesz(valuesz<minz)=minz;
    cmapz=mat2gray(valuesz)*128; cmapz=round(cmapz); cmapz(cmapz<1)=1; cmapz(cmapz>128)=128;
    lowcode=[0 0 1]; highcode=[1 1 0]; midcode=mean([lowcode;highcode]);
    cmapmid=round(128*(midz-minz)/(maxz-minz));
    cmap=makecmap_midpoint(lowcode,midcode,highcode,cmapmid);
    colorsz=cmap(cmapz,:);
    scatter(valuesx,valuesy,50,colorsz,'fill');
end
%titlestring='p21';
%axis([4 11 5.5 8.5]); %med
%axis([6 12 7 9.5]); %max
%axis([6 11 4.5 9]); %mean
%axis([7 9 4 7]);
axis([4.5 10 5 7.5]);
xstring='p21 log2(medRFU)'; ystring='CycD1 log2(medRFU)';
xlabel(xstring); ylabel(ystring);
set(gcf,'color','w','PaperPosition',[0 0 5 6]); %6 4.5
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
Hoechstval=IFdata(:,4);
%Ab1val=IFdata(:,5)-IFdata(:,14);
Ab1val=IFdata(:,8);
Ab1val(Ab1val<1)=1; Ab1val=log2(Ab1val);
%Ab2val=IFdata(:,6)-IFdata(:,15);
Ab2val=IFdata(:,9);
Ab2val(Ab2val<1)=1; Ab2val=log2(Ab2val);
Ab3val=log2(IFdata(:,13)); %max:13
Ab3val(Ab3val<1)=1; Ab3val=log2(Ab3val);
%%% Gate by Hoechst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hist(Hoechstval,0:10000:1000000);
G1cells=find(Hoechstval<250000); %pRbAbChar 20131217 Bgsub G1:<250000 G2:400000<x<600000
%G1cells=find(Hoechstval>400000 & Hoechstval<600000);
%G1cells=find(Hoechstval<225000 & Hoechstval>100000); %G1 default 200000
%G1cells=find(Hoechstval>300000); %G2
%G1cells=find(Hoechstval<600000); %All
Ab1=Ab1val(G1cells);
Ab2=Ab2val(G1cells);
Ab3=Ab3val(G1cells);
end
