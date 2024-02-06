function Stain_HistComparison(conditions,datadir)
EdUMeasured=0;
DisplayOption='EachVsCtrl'; %'Overlay' 'EachVsCtrl' '3D' 'Subtract'
showylabel=0;
%condnum=size(conditions,1);
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);
bmin=0; %5/6:0  5:2
bmax=4; %5/6:3  5:13
bstep=(bmax-bmin)/25; %pRb:50 pRb/tRb:100
bin=bmin:bstep:bmax;
[~,idx1]=min(abs(bin-1));
numbins=numel(bin);
binfill=[bin fliplr(bin)];
namecount=cell(uniquecondnum,1);
colorcode='krgcmb';
%colorcode=distributecolormap(jet,uniquecondnum); %alt:cool
%colors=colorcode+0.3*(1-colorcode);
colors=colorcode;
%scattercolor='br';
% c1=[1 0 0]; c2=[0 1 0]; c3=[0 0 1];
% color1=c1+.7*(1-c1); color2=c2+.7*(1-c2); color3=c3+.7*(1-c3);
% cmat={color1,color2,color3};
fontsizevar=8; %2col:8 4col:6
if showylabel
    mlvar=0.25; %default 0.25
    %pphvar=2; %4col per slide
    pphvar=3;
    %pphvar=4; %2col per slide
    %phvar=6; %1col per slide
else
    mlvar=0.1;
    %pphvar=1.75; %4col per slide
    pphvar=2.5;
    %pphvar=3.5; %2col per slide
end

pmat=cell(uniquecondnum,1);
for i=1:uniquecondnum
    Alldata=[];
    condrow=find(ismember(conditions(:,1),uniquenames{i}));
    for c=condrow'
        rowmat=cell2mat(conditions(c,2));
        colmat=cell2mat(conditions(c,3));
        sitemat=cell2mat(conditions(c,4));
        for row=rowmat
            for col=colmat
                for site=sitemat
                    %shot=wellnum2str(row,col,site);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    %load([datadir,shot,'.mat'],'IFdata');
                    load([datadir,shot,'.mat'],'IFdata');
                    Alldata=[Alldata;IFdata];
                end
            end
        end
    end
    Hoechstval=Alldata(:,3).*Alldata(:,4);
    c=1;
    %c=i;
    %Serum Release: 0nM,10nM,100nM,1uM
    G1minH=350000;
    G1maxH=550000;
    G2minH=850000;
    G2maxH=1100000;
    G1minE=12.5; G1maxE=15;
    G2minE=14; G2maxE=15.5;
    SminH=400000; SmaxH=1100000; SminE=18; SmaxE=22;
    %     %Cycling Cells: 0nM,10nM,100nM,1uM
%     %%%0nM,10nM,100nM,1uM
%     %G1minH=[250000,200000,200000,200000]; %1000000
%     %G1maxH=[450000,400000,400000,375000]; %2250000
%     %G2minH=[650000,550000,550000,500000];
%     %G2maxH=[950000,800000,800000,750000];
    xylim=[G1minH(c)-100000 G2maxH(c)+100000 G1minE-1 SmaxE+1];
    if EdUMeasured==1
        EdUval=log2(Alldata(:,3).*Alldata(:,9));
        Allcells=Hoechstval>xylim(1) & Hoechstval<xylim(2) & EdUval>xylim(3) & EdUval<xylim(4);
        G1cells=Hoechstval>G1minH(c) & Hoechstval<G1maxH(c) & EdUval>G1minE & EdUval<G1maxE;
        Scells=Hoechstval>SminH & Hoechstval<SmaxH & EdUval>SminE & EdUval<SmaxE;
        G2cells=Hoechstval>G2minH(c) & Hoechstval<G2maxH(c) & EdUval>G2minE & EdUval<G2maxE;       
    else
        Allcells=Hoechstval>0 & Hoechstval<800000;
        G1cells=Hoechstval>G1minH & Hoechstval<G1maxH;
        G2cells=Hoechstval>G1maxH;
    end
    Abvals=Alldata(:,4).*Alldata(:,3);
    %Abvals=Alldata(:,10)./Alldata(:,9);
    %Abvals=Alldata(:,8).*Alldata(:,3);
    Abvalsx=Alldata(:,5);
    %Abvalsx=Alldata(:,9)./Alldata(:,6);
    %xcoor=Alldata(:,1); ycoor=Alldata(:,2);
    posval=Abvals>=1; Abvals(Abvals<1)=1; Abvals=log2(Abvals);
    posvalx=Abvalsx>=1; Abvalsx(Abvalsx<1)=1; Abvalsx=log2(Abvalsx);
    
    %posval=Abvals>=1 & Alldata(:,6)>150; %100
    %Abvals(Abvals<1)=1; Abvals=log2(Abvals);
    
    %posval=Abvals>-1;
    %posval=Abvals>-1 & Alldata(:,5)>200;
    %%% compare subpopulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cellgating=Allcells;
    %histogramofsubsets(Abvals,{G1cells,Scells,SG2cells,G2cells},{'G1','S','SG2','G2'},bin);
    %histogramofsubsets(Abvals,{G1cells,Scells,G2cells},{'G1','S','G2'},bin);
    Abvals=Abvals(Cellgating & posval & posvalx);
    Abvalsx=Abvalsx(Cellgating & posval & posvalx);
    %EdUval=EdUval(Allcells & posval);
    %Abvalsx=Abvalsx(G1cells & posval);
    %xcoor=xcoor(G1cells & posval); ycoor=ycoor(G1cells & posval);
    %%% plot pdf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dscatter(Abvals,EdUval); %axis([0 13 6 22]); saveas(gcf,'h:\Downloads\Fig.jpg');
    %scatter(Abvalsx,Abvals,scattercolor(i),'o');
    %scatter(Abvals,Abvalsx,scattercolor(i),'o');
    %figure;
    %scatter(xcoor,ycoor,50,Abvals,'fill');
    %hist(Abvals,50); xlabel('pRb(p-S807/S811) log2(meanRFU)'); ylabel('# cells'); set(gcf,'color','w','PaperPosition',[0 0 4 4]); saveas(gcf,'h:\Downloads\Fig.jpg');
    %hist(Abvalsx,0:0.033:1.3); xlim([0 1.25]); xlabel('DHB (mean(cytoRFU)/mean(nucRFU))'); ylabel('# cells'); set(gcf,'color','w','PaperPosition',[0 0 4 4]); saveas(gcf,'h:\Downloads\Fig.jpg');

    %dscatter(Abvalsx,Abvals); axis([0.15 1.3 5 12]); xlabel('DHB ratio'); ylabel('pRb log2(medRFU)'); set(gcf,'color','w','PaperPosition',[0 0 4 4]); saveas(gcf,'h:\Downloads\Fig.jpg'); %All
    %dscatter(Abvalsx,Abvals); axis([0.15 1.5 5 12]); xlabel('DHB ratio'); ylabel('pRb log2(medRFU)'); set(gcf,'color','w','PaperPosition',[0 0 4 4]); saveas(gcf,'h:\Downloads\Fig.jpg'); %G2 0th min
    dscatter(Abvalsx,Abvals); xlabel('TxRed'); ylabel('EdU');
    pdfvals=histc(Abvals,bin);
    pdfvals=100*pdfvals/sum(pdfvals);
    %pdfvals=log2(pdfvals);pdfvals(pdfvals<0)=0;
    cdfvals=cumsum(pdfvals);
    cdfvals=100*cdfvals/max(cdfvals);
    cdfvals(cdfvals==max(cdfvals))=99;
    if strcmp(DisplayOption,'EachVsCtrl') && i==1
        pdfctrlvals=pdfvals;
        continue;
    end
    switch DisplayOption
        case 'Overlay'
            %fill(binfill,[pdfvals;zeros(numbins,1)],colors(i,:),'edgecolor',colors(i,:),'FaceAlpha', 0.4);
            %fill(binfill,[pdfvals;zeros(numbins,1)],[1 1 1],'edgecolor',colors(i,:));
            line(bin,pdfvals,'linewidth',2,'color',colors(i));
            
            %fill(binfill,[pdfvals;zeros(numbins,1)],'b','FaceAlpha', 0.3);
            %counterfill(bin,pdfvals,bstep,0.75,1.5,'r');
            %stairs(bin,c_elements,'color',colors(i,:),'linewidth',2);
            %ylim([0 5]);
            hold on;
        case 'EachVsCtrl'
            %subaxis(uniquecondnum-1,1,i-1,'ML',mlvar,'MR',0.02,'MT',0.02,'MB',0.1,'SH',0.02); %MB 0.15
            subaxis(uniquecondnum-1,1,i-1,'ML',mlvar,'MR',0.03,'MT',0.02,'MB',0.1,'SH',0.02); %MB 0.15
            hold on;
            %fill(binfill,[pdfctrlvals;zeros(numbins,1)],'b','FaceAlpha', 0.3);
            %fill(binfill,[pdfvals;zeros(numbins,1)],'r','FaceAlpha', 0.3);
            if i==2
                line(bin,pdfctrlvals,'linewidth',2,'color','k');
            end
            line(bin,pdfvals,'linewidth',2,'color','r');
            %fill(binfill,[pdfctrlvals;pdfvals(end:-1:1)],'r','facealpha',0.3); %difference
            %a=pdfvals-pdfctrlvals; [~,idxhigh]=min(abs(bin-1.5)); suma=sum(a(idxhigh:end));
            %counterfill_twocurves(bin,pdfvals,pdfctrlvals,bstep,1.5,3.0,'b');
            
            %meanvals=mean([pdfctrlvals pdfvals],2);
            %fill(binfill,[meanvals;zeros(numbins,1)],'b','FaceAlpha', 0.3); %mean
            %counterfill_zerobottom(bin,meanvals,bstep,0.75,1.5,'r');
            xlim([bmin bmax]);
            ylim([0 5]);
            if showylabel
                y=ylabel(char(uniquenames(i)),'rot',0,'fontsize',fontsizevar);
                set(y,'Units','Normalized','Position',[-0.2,0.4,0]);
            end
            100-cdfvals(idx1)
        case 'Scatter'
            subaxis(1,uniquecondnum,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02); %MB 0.15
            dscatter(Abvalsx,Abvals);
            %axis([0 1.2 5 12]); %SR
            %axis([0 1 5 13]); %CC
            hold on;
            %%% draw cell cycle gates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            rectangle('Position',[G1minH(c),G1minE,G1maxH(c)-G1minH(c),G1maxE-G1minE],'EdgeColor','r','linewidth',2);
            rectangle('Position',[SminH,SminE,SmaxH-SminH,SmaxE-SminE],'EdgeColor','r','linewidth',2);
            rectangle('Position',[G2minH(c),G2minE,G2maxH(c)-G2minH(c),G2maxE-G2minE],'EdgeColor','r','linewidth',2);
            axis(xylim);
            title(char(uniquenames(i)));
        case '3D'
            fill3(binfill',ones(size(binfill'))*-i,[pdfvals;zeros(numbins,1)],colors(i,:),'edgecolor',colors(i,:),'FaceAlpha',0.4);
            hold on;   
    end
    set(gca,'fontsize',fontsizevar);
    namecount{i}=char(uniquenames(i));
    pmat{i}=Abvals;
end
%p=ranksum(pmat{1},pmat{2}); fprintf('mean(set1)=%8.2f  mean(set2)=%8.2f\n',mean(pmat{1}),mean(pmat{2}));

%xlabel('pRb(807/811) log2(meanRFU)');
xlabel('pRb/tRb');
%xlabel('EdU log2(sumRFU)');
hxl=get(gca,'xlabel');
set(hxl,'fontsize',fontsizevar);

switch DisplayOption
    case 'Overlay'
        xlim([bmin bmax]);
        ylabel('pdf %');
        %ylabel('cdf');
        legend(char(namecount(:)),'location','northeast');
        set(gcf,'color','w','PaperPosition',[0 0 6 4]); %big:[8 6] med:[6 4] small:[3 2]
    case 'EachVsCtrl'
        set(gcf,'color','w','PaperPosition',[0 0 pphvar 4]);
    case 'Scatter'
        set(gcf,'color','w','PaperPosition',[0 0 12 3]);
    case '3D'
        xlim([bmin bmax]);
        zlabel('pdf (%)');
        set(gca,'ytick',[]);
        set(gca,'yticklabel',[]);
        view(25,45);
        %view(0,25);
        axis equal;
        axesLabelsAlign3D;
        set(gca,'position',[0.1 0.2 0.8 0.7]); %pRb
        %set(gca,'position',[-0.2 0.2 1 0.7]); %pRb/tRb
        set(gcf,'color','w','PaperPosition',[0 0 6 4]); %big:[8 6] med:[6 4] small:[3 2]
        set(gcf,'color','w','PaperPosition',[0 0 4 4]); %pRb:[4 4] pRb/tRb:[2 4]
        h_legend=legend(char(namecount(:)),'location','northeast');
        set(h_legend,'fontsize',6,'position',[0.75 0.75 0.1 0.1]); %pRb
        %set(h_legend,'fontsize',6,'position',[0.5 0.5 0.1 0.1]); %pRb/tRb
end

%textpos=[0.15 0.55 0.3 0.15];
%annotation('textbox',textpos,'String',['p-value = ',num2str(p)]);
saveas(gcf,'h:\Downloads\Fig.jpg');
%%% 
end