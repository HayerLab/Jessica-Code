function CellCount(conditions,rawpath,nucr,fold,frames)
%%% USAGE NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARG4: 0:no fold calculation 1:fold calculation
% ARG5: single frame (number) for no fold calc. 2-element vector [a b] for fold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucname='CFP_';
frame=frames(1);
if numel(frames)==2
    lastframe=frames(2);
end
condnum = size(conditions,1);
%%% get cell count per condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];
offset = [];
for i=1:condnum
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
                rawdir=[rawpath,shot,'\'];
                nuc_raw=log(single(imread([rawdir,nucname,num2str(frame),'.tif'])));
                nuc_mask=blobdetector(nuc_raw,nucr,-0.02); %default is -0.02.  Consider -0.03.
                [~,numcells(cc)]=bwlabel(nuc_mask);
                if fold
                    nuc_raw=log(single(imread([rawdir,nucname,num2str(lastframe),'.tif'])));
                    nuc_mask=blobdetector(nuc_raw,nucr,-0.02); %default is -0.02.  Consider -0.03.
                    [~,lastnum]=bwlabel(nuc_mask);
                    numcells(cc)=round(100*(lastnum-numcells(cc))/numcells(cc));
                end
            end
        end
    end
    data=[data;numcells];
    offset=[offset;i*ones(samplesize,1)];
end
keyboard;
clf;
boxplot(axes,data,offset,'labels',conditions(:,1),'labelorientation','inline');
ylabel('Cell Count');
title('Last Frame / First Frame','FontSize',16);
if fold
    ylabel('Fold Change (%)');
    title('Fold Change','FontSize',16);
end
set(gcf,'color','w','PaperPosition',[0 0 6 5]);
saveas(gcf,'h:\Downloads\Fig.jpg');