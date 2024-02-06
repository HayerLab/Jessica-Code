%function VerifyGating(row,col,site)
row='D';col='05';site='4';
%row='B';col='04';site='1';

projectpath='H:\Documents\Projects\';
experimentpath='2013-06-07_p21_cy2_deletions\Experiment_20130715\';
datadir=([projectpath,experimentpath,'Data\']);
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];
immunoframe=0;
%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([datadir,'tracedata_',shot,'.mat'],'tracedata','tracestats','jitters');
%%%%%% signal choice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal=2;
%%% calc stats %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigval=tracedata(:,:,signal);
sigval=sigval(:); sigval=sigval(~isnan(sigval));
sigdiff=diff(tracedata(:,:,signal),1,2);
abssigdiff=[sigdiff,zeros(size(sigdiff,1),1)];
relsigdiff=abssigdiff./tracedata(:,:,signal);
absdiff=abssigdiff(:); absdiff=absdiff(~isnan(absdiff));
reldiff=relsigdiff(:); reldiff=reldiff(~isnan(reldiff));
%%% plot statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samplesize=numel(absdiff);
cdfplot(sigval);
cdfplot(reldiff);
cdfplot(absdiff);
hist(reldiff,100);


%set(gcf,'color','w','PaperPosition',[0 0 15 10]); %3x4 or mini
%saveas(gcf,'h:\Downloads\Fig.jpg');
%close(gcf);