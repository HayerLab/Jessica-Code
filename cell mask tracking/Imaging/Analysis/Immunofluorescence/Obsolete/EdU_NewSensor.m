clear; close all;
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130719\';
datadir = [path,'Data\'];
%%%%%%%%%%%%%%%% NAME %%%%%%% ROWS %%%%%% COLS %%%% SITES %%%
conditions = {'DMSO 0.1%',   [2 3]  ,        7   ,  1:4};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'Units','pixels');                %sets screensize units by pixels
screendims = get(0,'ScreenSize');       %get screensize in pixels
screenx = screendims(3);
screeny = screendims(4);

condnum = size(conditions,1);
framesperhr = 5;
bstep = 1/framesperhr;
bmin = 0;
bmax = 15;
bin = [bmin:bstep:bmax];

stairfig = figure;
set(stairfig,'Position',[round(0.1*screenx) round(0.1*screeny) round(0.8*screenx) round(0.8*screeny)]); %3x3
set(stairfig,'color','w','PaperPosition',[0 0 12 8]); %cumulative stairs (2x3)
palette = 'krgc';


hold on;
row = cell2mat(conditions(1,2));
col = cell2mat(conditions(1,3));
site = cell2mat(conditions(1,3));
cd ..\Analysis
[sensor,EdU]=getEdUdata(row,col,site);

%%% visualize data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
cd ..\Functions
if correlationview
    b1min = 0; b1max = 15;
    b2min = 0; b2max = 15;
    subaxis(2,3,i,'ML',0.1,'MR',0.05,'MT',0.05,'MB',0.1,'SV',0.1);
    plot(values2,values,'.','markersize',8);
    %dscatter(values2,values,'PLOTTYPE','contour');
    title(char(conditions(i,1)));
    axis([b2min b2max b1min b1max]);
end
saveas(gcf,'h:\Downloads\Fig.jpg');
%close(histfig);
cd ..\Analysis; %return to this directory