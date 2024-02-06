function Hoechst_EdU_Ab(conditions,datadir)
condnum=size(conditions,1);
boxplotdata=[];
offset=[];
%%% Mouse log2
bmin=7;  %G1: med 7  max 7  int 16     S: med:7 max 7 sum 16    G2: med 7  max 7 sum 16
bmax=11; %G1: med 10 max 12 int 19     S: med:11 max 14 sum 20   G2: med 11 max 14 sum 21
%%% Goat raw
%bmin=7;    %G1: med 50     S: med 0    G2: med 0 max 100
%bmax=10;  %G1: med 350    S: med 400  G2: med 400 max 800
bstep=(bmax-bmin)/100; %pdf 30
bin=bmin:bstep:bmax;
namecount=cell(condnum,1);
colorcode='krgcmb';
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
    values=AbG1vals;
    %%% store boxplot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    boxplotdata=[boxplotdata;values];
    offset=[offset;i*ones(length(values),1)];
    %%% plot cdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_elements = histc(values,bin);
    c_elements = cumsum(n_elements);
    c_elements = 100*c_elements/max(c_elements);
    c_elements(c_elements==max(c_elements))=99;
    %n_elements= 100*n_elements/sum(n_elements);
    stairs(bin,c_elements,colorcode(i),'linewidth',2);
    namecount{i}=[char(conditions(i,1))];
end
legend(char(namecount(:)),'location','southeast');
xlabel('pRb(p-S780) medianRFU'); ylabel('cdf (%)'); %ylabel('count');
set(gcf,'color','w','PaperPosition',[0 0 6 4.5]); %6 4.5
saveas(gcf,'h:\Downloads\Fig.jpg');
%%% 
end


function [AbG1,AbS,AbG2]=main(datadir,shot)
load([datadir,shot,'.mat'],'IFdata');
Hoechstval=IFdata(:,4);
EdUval=log2(IFdata(:,10));
Abval=log2(IFdata(:,5));
%Abval=log2(IFdata(:,7).*IFdata(:,3));
%%% Gate by cell cycle phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Scells=find(Hoechstval<600000 & EdUval>8);  %600000 8
G1cells=find(Hoechstval<450000 & EdUval<7.5); %325000 7.5
G2cells=find(Hoechstval>400000 & Hoechstval<600000 & EdUval<7.5); %425000 600000 7
AbG1=Abval(G1cells); AbS=Abval(Scells); AbG2=Abval(G2cells);
end
