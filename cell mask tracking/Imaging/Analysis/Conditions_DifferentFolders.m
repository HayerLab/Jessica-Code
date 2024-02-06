%%% DEFINE CONDITIONS PER WELL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 2014-02-08_H2B-DHB-p21dCy1dK
%4:mean(Hoechst) 5:med(H2B) 6:med(DHBnuc) 7:med(p21dCy1dK) 8:mean(top50%(DHBring))
rowvar=4;
conditions={
    'Media',4,8,1,'H:\Documents\Projects\2014-02-08_H2B-DHB-p21dCy1dK\Data\';
    'Media',5,8,1,'H:\Documents\Projects\2014-02-08_H2B-DHB-p21dCy1dK\Data\';
};

%%% CHOOSE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% validated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stain_TrendComparison(conditions,datadir);
%Stain_HistComparison(conditions,datadir);
%Stain_3Dmap(conditions,datadir);
%Geminin(conditions,datadir);
%pipdeg_CDK2(conditions,datadir);
%Timelapse_IF_Cdk2(conditions,datadir);
%Visualize_CellCycleGating(conditions,datadir);
%Timelapse_General_CombineData(conditions,datadir);
Timelapse_General_CombineDataDifferentFolders(conditions);