function Stain_HistComparison_Interwell(conditions,datadir)
EdUMeasured=0;
condnum=size(conditions,1);
boxplotdata=[];
offset=[];
bmin=2; %780:7
bmax=8; %780:12
bstep=(bmax-bmin)/30; %pdf 30
bin=bmin:bstep:bmax;
binfill=[bin fliplr(bin)];
namecount=cell(condnum,1);
%colorcode='krgcmb';
colorcode=distributecolormap(jet,condnum); %alt:cool
colors=colorcode+0.3*(1-colorcode);

% c1=[1 0 0]; c2=[0 1 0]; c3=[0 0 1];
% color1=c1+.7*(1-c1); color2=c2+.7*(1-c2); color3=c3+.7*(1-c3);
% cmat={color1,color2,color3};

pmat=cell(condnum,1);
figure;hold on;
for i=1:condnum
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    Alldata=[];
    %Abvals=[];
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                %shot=wellnum2str(row,col,site);
                shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                load([datadir,shot,'.mat'],'IFdata');
                Alldata=[Alldata;IFdata];
            end
        end
    end
    Hoechstval=Alldata(:,3).*Alldata(:,4);
    Allcells=Hoechstval>150000 & Hoechstval<250000;
    if EdUMeasured==1
        EdUval=log2(Alldata(:,3).*Alldata(:,7));
        G1cells=Hoechstval>250000 & Hoechstval<550000 & EdUval> 13 & EdUval<16.5;
        Scells=Hoechstval>300000 & Hoechstval<700000 & EdUval>19 & EdUval<22;
        SG2cells=Hoechstval>600000 & Hoechstval<1000000 & EdUval>17.5 & EdUval<18.5;
        G2cells=Hoechstval>700000 & Hoechstval<1000000 & EdUval>14.5 & EdUval<17;
    else
        G1cells=find(Hoechstval<415000); %20131126:<275000
        G2cells=find(Hoechstval>415000);
    end
    Abvals=Alldata(:,5); Abvals(Abvals<1)=1; Abvals=log2(Abvals);
    %%% compare subpopulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %histogramofsubsets(Abvals,{G1cells,Scells,SG2cells,G2cells},{'G1','S','SG2','G2'},bin);
    %histogramofsubsets(Abvals,{G1cells,Scells,G2cells},{'G1','S','G2'},bin);
    Abvals=Abvals(Allcells);
    %%% plot pdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n_elements = histc(Abvals,bin);
    n_elements= 100*n_elements/sum(n_elements);
    c_elements = cumsum(n_elements);
    c_elements = 100*c_elements/max(c_elements);
    c_elements(c_elements==max(c_elements))=99;
    fill(binfill,[n_elements;zeros(length(n_elements),1)],colors(i,:),'edgecolor',colors(i,:),'FaceAlpha', 0.4);
    %stairs(bin,c_elements,colorcode(i),'linewidth',2);
    namecount{i}=char(conditions(i,1));
    pmat{i}=Abvals;
end
%p=ranksum(pmat{1},pmat{2}); fprintf('mean(set1)=%8.2f  mean(set2)=%8.2f\n',mean(pmat{1}),mean(pmat{2}));
legend(char(namecount(:)),'location','northwest');
%xlabel('pRb(807/811) log2(meanRFU)');
%xlabel('p21 log2(meanRFU)');
%xlabel('pRb(780) log2(meanRFU)');
xlabel('CycA log2(meanRFU)');
%xlabel('CycD1 log2(meanRFU)');
ylabel('pdf (%)'); %ylabel('count');
xlim([bmin bmax]);
%textpos=[0.15 0.55 0.3 0.15];
%annotation('textbox',textpos,'String',['p-value = ',num2str(p)]);
set(gcf,'color','w','PaperPosition',[0 0 6 4]); %big:[8 6] med:[6 4] small:[3 2]
saveas(gcf,'h:\Downloads\Fig.jpg');
saveas(gcf,'h:\Downloads\Fig.eps');
%%% 
end