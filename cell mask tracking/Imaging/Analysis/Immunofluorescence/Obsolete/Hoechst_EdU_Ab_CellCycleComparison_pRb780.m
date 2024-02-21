function Hoechst_EdU_Ab(conditions,datadir)
condnum=size(conditions,1);
boxplotdata=[];
offset=[];
bmin=7;  %med 7  max 7  int 15
bmax=10; %med 11 max 14 int 21
bstep=(bmax-bmin)/50;
bin=bmin:bstep:bmax;
namecount=cell(3,1);
colorcode='krgc';
figure;hold on;
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    AbG1vals=[];
    AbSvals=[];
    AbG2vals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [AbG1val,AbSval,AbG2val]=main(datadir,shot);
                AbG1vals=[AbG1vals;AbG1val];
                AbSvals=[AbSvals;AbSval];
                AbG2vals=[AbG2vals;AbG2val];
            end
        end
    end
end
values=AbG1vals;
boxplotdata=[boxplotdata;values];
offset=[offset;i*ones(length(values),1)];
n_elements = histc(values,bin);
n_elements= 100*n_elements/sum(n_elements);
stairs(bin,n_elements,colorcode(1),'linewidth',2);

namecount{1}='G1';
values=AbSvals;
boxplotdata=[boxplotdata;values];
offset=[offset;i*ones(length(values),1)];
n_elements = histc(values,bin);
n_elements= 100*n_elements/sum(n_elements);
stairs(bin,n_elements,colorcode(2),'linewidth',2);
namecount{2}='S';

values=AbG2vals;
boxplotdata=[boxplotdata;values];
offset=[offset;i*ones(length(values),1)];
n_elements = histc(values,bin);
n_elements= 100*n_elements/sum(n_elements);
stairs(bin,n_elements,colorcode(3),'linewidth',2);
namecount{3}='G2';

xlim([bmin bmax]);
legend(char(namecount(:)),'location','northeast');
%title('Cell Cycle comparison');
xlabel('pRb(p-S780) log2(maxRFU)'); ylabel('cdf (%)'); ylabel('pdf (%)')
set(gcf,'color','w','PaperPosition',[0 0 6 4.5]);
saveas(gcf,'h:\Downloads\Fig_cdfs.jpg');
%{
figure,boxplot(axes,boxplotdata,offset,'labels',namecount);
title('Cell Cycle comparison');ylabel('pRb(p-S780) ');
set(gcf,'color','w','PaperPosition',[0 0 4 5]);
saveas(gcf,'h:\Downloads\Fig_boxplots.jpg');
%} 
end


function [AbG1,AbS,AbG2]=main(datadir,shot)
load([datadir,shot,'.mat'],'IFdata');
Hoechstval=IFdata(:,4);
EdUval=log2(IFdata(:,10));
Abval=log2(IFdata(:,9));
%Abval=log2(IFdata(:,7).*IFdata(:,3));
%%% Gate by cell cycle phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Scells=find(Hoechstval<600000 & EdUval>9);
G1cells=find(Hoechstval<300000 & EdUval<7);
G2cells=find(Hoechstval>425000 & Hoechstval<600000 & EdUval<7);
AbG1=Abval(G1cells); AbS=Abval(Scells); AbG2=Abval(G2cells);
end