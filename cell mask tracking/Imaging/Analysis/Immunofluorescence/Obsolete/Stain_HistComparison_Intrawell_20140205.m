function Hoechst_EdU_Ab(conditions,datadir)
condnum=size(conditions,1);
boxplotdata=[];
offset=[];
bmin=2;  %med 2  max 7  int 15
bmax=12; %med 11 max 14 int 21
bstep=(bmax-bmin)/50;
bin=bmin:bstep:bmax;
binfill=[bin fliplr(bin)];
c1=[1 0 0]; c2=[0 1 0]; c3=[0 0 1];
color1=c1+.7*(1-c1); color2=c2+.7*(1-c2); color3=c3+.7*(1-c3);
cmat={color1,color2,color3};
namecount=cell(3,1);
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
%fill(binfill,[n_elements;ones(length(n_elements),1)*-1],cmat{1},'edgecolor',cmat{1},'FaceAlpha', 0.4);
fill(binfill,[n_elements;zeros(length(n_elements),1)],cmat{1},'edgecolor',cmat{1},'FaceAlpha', 0.4);

namecount{1}='G1';
values=AbSvals;
boxplotdata=[boxplotdata;values];
offset=[offset;i*ones(length(values),1)];
n_elements = histc(values,bin);
n_elements= 100*n_elements/sum(n_elements);
fill(binfill,[n_elements;zeros(length(n_elements),1)],cmat{2},'edgecolor',cmat{2},'FaceAlpha', 0.4);
namecount{2}='S';

values=AbG2vals;
boxplotdata=[boxplotdata;values];
offset=[offset;i*ones(length(values),1)];
n_elements = histc(values,bin);
n_elements= 100*n_elements/sum(n_elements);
fill(binfill,[n_elements;zeros(length(n_elements),1)],cmat{3},'edgecolor',cmat{3},'FaceAlpha', 0.4);
namecount{3}='G2';

xlim([bmin bmax]);
legend(char(namecount(:)),'location','northeast');
%title('Cell Cycle comparison');
%xlabel('pAKT(S473) log2(medianRFU)'); ylabel('pdf (%)');
xlabel('pERK(T202/Y204) log2(medianRFU)'); ylabel('pdf (%)');
set(gcf,'color','w','PaperPosition',[0 0 8 6]);
saveas(gcf,'h:\Downloads\Fig.jpg');
saveas(gcf,'h:\Downloads\Fig.eps');
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
Abval=log2(IFdata(:,5));
%Abval=log2(IFdata(:,7).*IFdata(:,3));
%%% Gate by cell cycle phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1cells=find(Hoechstval>200000 & Hoechstval<500000 & EdUval>7.5 & EdUval<8.5);
Scells=find(Hoechstval>300000 & Hoechstval<900000 & EdUval>10 & EdUval<12.5);
G2cells=find(Hoechstval>620000 & Hoechstval<1000000 & EdUval>7.8 & EdUval<8.8);
AbG1=Abval(G1cells); AbS=Abval(Scells); AbG2=Abval(G2cells);
end
