function ShadingEval(conditions,datadir)
EdUMeasured=0;
%condnum=size(conditions,1);
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);
bmin=14; %5/6:0     5:4  sum5:13  6:6  sum6:14
bmax=21; %5/6:2.25  5:13  sum5:21  6:12 sum6:21
bstep=(bmax-bmin)/40; %pdf 30
bin=bmin:bstep:bmax;
binfill=[bin fliplr(bin)];
namecount=cell(uniquecondnum,1);
%colorcode='krgcmb';
colorcode=distributecolormap(jet,uniquecondnum); %alt:cool
colors=colorcode+0.3*(1-colorcode);
%scattercolor='br';
% c1=[1 0 0]; c2=[0 1 0]; c3=[0 0 1];
% color1=c1+.7*(1-c1); color2=c2+.7*(1-c2); color3=c3+.7*(1-c3);
% cmat={color1,color2,color3};

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
                    load([datadir,shot,'.mat'],'IFdata');
                    Alldata=[Alldata;IFdata];
                end
            end
        end
    end
    Hoechstval=Alldata(:,3).*Alldata(:,4);
    subprctile=25; %lower/upper percentile to compare
    [peakstart,lowerthresh,peakmid,upperthresh,peakend,lastval]=ParameterizeFirstPeak(Hoechstval,subprctile);
    if EdUMeasured==0
        G1cells=Hoechstval>peakstart & Hoechstval<peakend; %[350000 550000]
    else
        EdUval=log2(Alldata(:,3).*Alldata(:,7));
        G1cells=Hoechstval>peakstart & Hoechstval<peakend & EdUval> 11 & EdUval<14.5;
    end
    xcoor=Alldata(:,1); ycoor=Alldata(:,2);
    highvals=Hoechstval>upperthresh; %Corrected:375000  NotCorrected:475000
    lowvals=Hoechstval<lowerthresh;  %Corrected:325000  NotCorrected:400000
    %hint=Alldata(:,4);
    %highvals=hint>1100; %SC:800 SNC:1100 SNC(area):600
    %lowvals=hint<700; %SC:600 SNC:700 SNC(area):400
    
    highcells=G1cells & highvals; lowcells=G1cells & lowvals;
    Hoechsthigh=Hoechstval(highcells);
    Hoechstlow=Hoechstval(lowcells);
    %Hoechsthigh=hint(highcells);
    %Hoechstlow=hint(lowcells);
    xcoorhigh=xcoor(highcells); ycoorhigh=ycoor(highcells);
    xcoorlow=xcoor(lowcells); ycoorlow=ycoor(lowcells);
    %%% Aspect Comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    totalwidth=lastval-peakstart;
    peakwidth=peakend-peakstart;
    relpeakwidth=100*peakwidth/totalwidth;
    fprintf('peakwidth = %0.0f\n',peakwidth);
    fprintf('relpeakwidth = %0.0f%%\n',relpeakwidth);
    %%% spatial comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure,scatter(xcoorhigh,ycoorhigh,50,Hoechsthigh,'fill');
    figure,scatter(xcoorlow,ycoorlow,50,Hoechstlow,'fill');
    
    %namecount{i}=char(uniquenames(i));
end
%legend(char(namecount(:)),'location','northeastoutside');
set(gcf,'color','w','PaperPosition',[0 0 4 4]); %big:[8 6] med:[6 4] small:[3 2]
saveas(gcf,'h:\Downloads\Fig.jpg');
end