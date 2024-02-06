clear; close all;
codepath='h:\Documents\Timelapse\Code\Development\';
cd([codepath,'Functions\']); %change directory for function calls
path = 'H:\Documents\Experiments\2013-06-07_p21_cy2_deletions\Imaging_20130719\';
datadir = [path,'Data\'];
shot = 'B_07_1';
load([datadir,'wellsss_IF_',shot,'.mat'],'wellsss');
stainframe=size(wellsss,3);
Cy5data=wellsss{stainframe}(:,8);
RFPdata=wellsss{stainframe}(:,6);
%plot(Cy5data,RFPdata,'.','markersize',8);
dscatter(Cy5data,RFPdata,'PLOTTYPE','contour');
title('EdU vs d57 sensor');
axis([0 2.5 0 3]);
keyboard;
%saveas(gcf,'h:\Downloads\Fig.jpg');
%close(histfig);
cd ..\Analysis; %return to this directory