function histogramofsubsets(Abvals,subsets,names,bin)
figure; hold on;
numsubsets=numel(subsets);
colorcode=distributecolormap(jet,numsubsets); %alt:cool
colors=colorcode+0.3*(1-colorcode);
binfill=[bin fliplr(bin)];
for i=1:numsubsets
    values=Abvals(subsets{i});
    n_elements = histc(values,bin);
    n_elements= 100*n_elements/sum(n_elements);
    fill(binfill,[n_elements;zeros(length(n_elements),1)],colors(i,:),'edgecolor',colors(i,:),'FaceAlpha', 0.4);
end

xlim([min(bin) max(bin)]);
legend(char(names(:)),'location','northeast');
xlabel('CycA log2(meanRFU)');
%xlabel('CycD1 log2(meanRFU)');
ylabel('pdf (%)');
set(gcf,'color','w','PaperPosition',[0 0 8 6]);
saveas(gcf,'h:\Downloads\FigIntrawell.jpg');