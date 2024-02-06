imagepath='H:\Images\';
experimentpath='2013-08-05_p21_cy1_deletions\Experiment_20130831\';
nucr=8;
nucname='CFP_';
firstframe=1;
lastframe=234;
%%% define experimental wells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conditions={
    'No Sensor',2:4,2,1:2;
    'WT',2:4,3,1:2;
    'N',2:4,4,1:2;
    'mCDK',2:4,5,1:2;
    'mPCNA',2:4,6,1:2;
    'mCDKmPCNA',2:4,7,1:2;
    'mCy2',2:4,8,1:2;
    'C',2:4,9,1:2;
    'C-mCy2',2:4,10,1:2;
};
condnum = size(conditions,1);
%%% get cell count per condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];
offset = [];
for i=1:condnum
    hold on;
    rowmat=cell2mat(conditions(i,2));
    colmat=cell2mat(conditions(i,3));
    sitemat=cell2mat(conditions(i,4));
    samplesize=numel(rowmat)*numel(colmat)*numel(sitemat);
    numcells=zeros(samplesize,1);
    cc=0;
    for row=rowmat
        for col=colmat
            for site=sitemat
                cc=cc+1;
                shot=wellnum2str(row,col,site);
                rawdir=[imagepath,experimentpath,'Raw\',shot,'\'];
                nuc_raw=log(single(imread([rawdir,nucname,num2str(firstframe),'.tif'])));
                nuc_mask=blobdetector(nuc_raw,nucr,-0.03); %default is -0.02.  Consider -0.03.
                [~,firstnum]=bwlabel(nuc_mask);
                nuc_raw=log(single(imread([rawdir,nucname,num2str(lastframe),'.tif'])));
                nuc_mask=blobdetector(nuc_raw,nucr,-0.03); %default is -0.02.  Consider -0.03.
                [~,lastnum]=bwlabel(nuc_mask);
                numcells(cc)=round(100*(lastnum-firstnum)/firstnum);
            end
        end
    end
    data=[data;numcells];
    offset=[offset;i*ones(samplesize,1)];
end
clf;
boxplot(axes,data,offset,'labels',conditions(:,1),'labelorientation','inline');
ylabel('Fold Change (%)'); title('Last Frame / First Frame','FontSize',16);
set(gcf,'color','w','PaperPosition',[0 0 6 5]);
saveas(gcf,'h:\Downloads\Fig.jpg');