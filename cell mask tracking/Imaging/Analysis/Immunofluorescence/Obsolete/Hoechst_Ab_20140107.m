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
bmin=0; %6 CycD 4
bmax=9; %10 CycD 11
bstep=(bmax-bmin)/50; %pdf 30
bin=bmin:bstep:bmax;
namecount=cell(condnum,1);
sigmean=ones(1,condnum)*NaN;
sigstd=ones(1,condnum)*NaN;
sigtime=0:11;
%colorcode='krgcmb'; colorcode='krc';
colorcode=[0 0 0; 1 0 0; 0 1 0; 0 0 1];
%colorcode=[0 0 0; 0 0 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1];
%colorcode=[0 0 0; 0 0 0; 1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1];
%colorcode=distributecolormap(jet,condnum);
figure;hold on;
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
    values1=Ab2vals;
    %%% store boxplot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %boxplotdata=[boxplotdata;values1];
    %offset=[offset;i*ones(length(values1),1)];
    %%% plot cdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_elements = histc(values1,bin);
    c_elements = cumsum(n_elements);
    c_elements = 100*c_elements/max(c_elements);
    c_elements(c_elements==max(c_elements))=99;
    n_elements= 100*n_elements/sum(n_elements);
    stairs(bin,c_elements,'color',colorcode(i,:),'linewidth',2);
    namecount{i}=[char(conditions(i,1))];
    %[sigmean(i),sigstd(i)]=normfit(values1);
    sigmean(i)=nanmean(values1);
    sigstd(i)=nanstd(values1);
end
titlestring='p21';
descstring='p21 log2(meanRFU)';
legend(char(namecount(:)),'location','northwest');
xlabel(descstring); ylabel('cdf (%)'); %ylabel('count');
set(gcf,'color','w','PaperPosition',[0 0 9 6]); %6 4.5
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
%Ab1val=log2(IFdata(:,8).*IFdata(:,3)); %Ab1val=log2(IFdata(:,5));
Ab1val=log2(IFdata(:,8));
Ab2val=log2(IFdata(:,9));
Ab2val(Ab2val<1)=1; Ab2val=log2(Ab2val);
Ab3val=log2(IFdata(:,10));
%%% Gate by Hoechst %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1cells=find(Hoechstval<225000); %cycDp21SerumRelease: 300000
%G2cells=find(Hoechstval>350000 & Hoechstval<550000);
Ab1=Ab1val(G1cells);
Ab2=Ab2val(G1cells);
Ab3=Ab3val(G1cells);
end
