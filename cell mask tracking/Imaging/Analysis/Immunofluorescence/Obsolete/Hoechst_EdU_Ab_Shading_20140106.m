function Hoechst_EdU_Ab(conditions,datadir)
EdUMeasured=0;
condnum=size(conditions,1);
boxplotdata=[];
offset=[];
%%% pAKT
bmin=4;  %G1: med 3  max 7  int 16     S: med:7 max 7 sum 16    G2: med 7  max 7 sum 16
bmax=7; %G1: med 9 max 12 int 19     S: med:11 max 14 sum 20   G2: med 11 max 14 sum 21
%%% pERK
%bmin=7;    %G1: med 50     S: med 0    G2: med 0 max 100
%bmax=10;  %G1: med 350    S: med 400  G2: med 400 max 800
bstep=(bmax-bmin)/100; %pdf 30
bin=bmin:bstep:bmax;
namecount=cell(condnum,1);
colorcode='krgcmb';
%colorcode=distributecolormap(cool,condnum);
binfill=[bin fliplr(bin)];
c1=[1 0 0]; c2=[0 0 1]; c3=[0 0 1];
color1=c1+.7*(1-c1); color2=c2+.7*(1-c2); color3=c3+.7*(1-c3);
cmat={color1,color2,color3};
figure;hold on;
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    Ab1Allvals=[];
    Ab1G1vals=[];
    Ab1Svals=[];
    Ab1G2vals=[];
    Ab2Allvals=[];
    Ab2G1vals=[];
    Ab2Svals=[];
    Ab2G2vals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                [Ab1Allval,Ab1G1val,Ab1Sval,Ab1G2val,Ab2Allval,Ab2G1val,Ab2Sval,Ab2G2val]=main(datadir,shot,EdUMeasured);
                Ab1Allvals=[Ab1Allvals;Ab1Allval];
                Ab1G1vals=[Ab1G1vals;Ab1G1val];
                Ab1Svals=[Ab1Svals;Ab1Sval];
                Ab1G2vals=[Ab1G2vals;Ab1G2val];
                Ab2Allvals=[Ab2Allvals;Ab2Allval];
                Ab2G1vals=[Ab2G1vals;Ab2G1val];
                Ab2Svals=[Ab2Svals;Ab2Sval];
                Ab2G2vals=[Ab2G2vals;Ab2G2val];
            end
        end
    end
    values1=Ab1G1vals;
    values2=Ab2G1vals;
    %%% store boxplot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %boxplotdata=[boxplotdata;values1];
    %offset=[offset;i*ones(length(values1),1)];
    %%% plot cdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_elements = histc(values1,bin);
    c_elements = cumsum(n_elements);
    c_elements = 100*c_elements/max(c_elements);
    c_elements(c_elements==max(c_elements))=99;
    n_elements= 100*n_elements/sum(n_elements);
    %stairs(bin,n_elements,colorcode(i),'linewidth',2);
    fill(binfill,[n_elements;zeros(length(n_elements),1)],cmat{i},'edgecolor',cmat{i},'FaceAlpha', 0.4);
    namecount{i}=[char(conditions(i,1))];
%     if i==1
%         a=values1;
%     else
%         b=values1;
%     end
end
%p=ranksum(a,b);
legend(char(namecount(:)),'location','northwest');
%xlabel('pAKT(S473) log2(medianRFU)');
%xlabel('pERK(T202/Y204) log2(medianRFU)');
xlabel('pRb(807/811) log2(meanRFU)');
ylabel('pdf (%)'); %ylabel('count');
%textpos=[0.15 0.55 0.3 0.15];
%annotation('textbox',textpos,'String',['p-value = ',num2str(p)]);
set(gcf,'color','w','PaperPosition',[0 0 4 3]); %6 4.5
saveas(gcf,'h:\Downloads\Fig.jpg');
saveas(gcf,'h:\Downloads\Fig.eps');
%%% 
end


function [Ab1All,Ab1G1,Ab1S,Ab1G2,Ab2All,Ab2G1,Ab2S,Ab2G2]=main(datadir,shot,EdUMeasured)
load([datadir,shot,'.mat'],'IFdata');
Hoechstval=IFdata(:,4);
if EdUMeasured==1
    EdUval=log2(IFdata(:,13));
end
Ab1val=log2(IFdata(:,8));
%Abval=log2(IFdata(:,7).*IFdata(:,3));
Ab2val=log2(IFdata(:,6));
%%% Gate by cell cycle phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Allcells=find(Hoechstval>200000 & Hoechstval<1100000);
if EdUMeasured==1
    G1cells=find(Hoechstval>200000 & Hoechstval<500000 & EdUval>7.5 & EdUval<8.5);
    Scells=find(Hoechstval>300000 & Hoechstval<900000 & EdUval>10 & EdUval<12.5);
    G2cells=find(Hoechstval>620000 & Hoechstval<1000000 & EdUval>7.8 & EdUval<8.8);
else
    G1cells=find(Hoechstval<275000);
    Scells=G1cells;
    G2cells=find(Hoechstval>275000);
end
Ab1All=Ab1val(Allcells); Ab1G1=Ab1val(G1cells); Ab1S=Ab1val(Scells); Ab1G2=Ab1val(G2cells);
Ab2All=Ab2val(Allcells); Ab2G1=Ab2val(G1cells); Ab2S=Ab2val(Scells); Ab2G2=Ab2val(G2cells);
end
