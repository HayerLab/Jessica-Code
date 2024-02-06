%%% DEFINE CONDITIONS PER WELL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 20131126 R-point SR\20140317 3C CDK4i Wee1i
% conditions={
%     'CDK4/6i',5:7,5:6,1:2;
% };

%%% 20131209 Pipdeg Optimization\20140402 PFACS
% conditions={
%     'Ctrl',1:8,10,1:4;
% };

%%% 20131126 R-point SR\20140317 3C CDK4i Wee1i
% conditions={
%     '2i',5,3,1;
%     };

%%% 20131213 R-point CC\20140318 CyclingCDK4iWee1iCDC25i
% conditions={
%     'Ctrl',2,4,1;
% };

%%% Heewon
% conditions={
%     'Ctrl',4,3,1;
% };

%%% 20131213 R-point CC\20140405 2C siRNA CDK4i
% conditions={
%     %'Ctrl',1:8,11,1:4;
%     %'CDK4/6i',1:8,12,1:4;
%     %'siCtrl',1:8,1:2,1:4;
%     %'siCtrl+4/6i',1:8,3:4,1:4;
%     %'CycE',1:8,5:6,1:4;
%     'siCycD1',1:8,10,1:4;
%     };

%%% 20131126 R-point SR\20140404 siCycE siCycD CDK4i CDK2i
% cv=0;
% conditions={
%     %'Ctrl',1:4,5+cv,1;
%     %4/6iconc + 2i60uM
%     %'1nM',1:4,5+cv,1;
%     %'10nM',1:4,6+cv,1;
%     %'100nM',1:4,7+cv,1;
%     %'1uM',1:4,8+cv,1;
%     
%     %'DMSO 9hrs',1:4,3,1;
%     'siCtrl 10hrs',5:8,1,1;
%     'siCtrl 4/6i',5:8,2,1;
%     %'siCtrl 2i',5:8,3,1;
%     %'siCtrl 4/6i+2i',5:8,4,1;
%     %'siCycE',5:8,5,1;
%     %'siE2F1',5:8,7,1;
%     %'siCycD1',5:8,11,1;
%     };

%%% 20131213 R-point CC\20140407 2C CDK2 Hysteresis
% conditions={
%     %'Ctrl',1:8,4,1:4;
%     '4/6i',1:8,2,1:4;
%     };

%%% 20131213 R-point CC\20140406 siRNA CDK4i
% conditions={
%     'DMSO 0.3%',5:8,5,1;
%     %'4/6i',5:8,6,1;
%     %'2i',5:8,7,1;
%     %'4/6i + 2i',5:8,8,1;
%     
%     %'siCycE',5:8,9,1;
%     %'siCycE+4/6i',5:8,10,1;
%     'siCDK4/6',5:8,12,1;
%     };

%%% 20131213 R-point CC\20140413 2C DHB Hysteresis
% conditions={
%     %'2i wash2hr',2:7,3,1:8;
%     %'2i+4/6i wash2hr',2:7,4,1:8;
%     %'2i wash3hr',2:7,5,1:8;
%     %'2i+4/6i wash3hr',2:7,6,1:8;
%     
%     'DMSO',2:7,7,1:8;
%     %'4/6i',2:7,8,1:8;
%     %'2i',2:7,9,1:8;
%     %'4/6i + 2i',2:7,10,1:8;
%     };

%%% 20131126 R-point SR\20140415 CDK46hysteresis-under-CDK2i
% rv=[9]; sv=1:9;
% conditions={
%     %%% EdU
%     %'17hr',12:13,1,sv;
%     %'11.5hr',12,2:3,sv;
%     
%     %%% pRb/tRb
%     %'DMSO 0.1%',rv,2,sv;
%     %'4/6i 1uM',rv,3,sv;
%     %'4/6i + 2i',3:4,10,sv;
%     %'4/6i + siCycE',7:8,10,sv;
%     %'4/6i + siE2F1',5:8,11,sv;
%     
%     %%% 4/6Hys-under-2i
%     'Ctrl',rv,4,sv;
%     '0nM',rv,4,sv;
%     '3nM',rv,5,sv;
%     '10nM',rv,6,sv;
%     '30nM',rv,7,sv;
%     '100nM',rv,8,sv;
%     '300nM',rv,9,sv;
%     '1uM',rv,10,sv;
%     };

%%% 20131213 R-point CC\20131217 3Doldrabbit XPcdkpanel
% conditions={
%     'Ctrl',2:4,1,1:4;
%     };

%%% Heewon breakbeforemitosis
% conditions={
%     'ctrl',2,4,1;
%     };

%%% pipdeg pseudoFACS retry
% conditions={
%     'ctrl',7,2:6,1:4;
%     };

%%% 20131213 R-point CC\20140423 3C G1S 46i
% conditions={
%     'sample',2,5,2;
%     %'DMSO 0.1%',2,2:6,1:4;
%     %'4/6i',3,2:6,1:4;
%     %'Cdc25i 10uM',4,2:6,1:4;
%     %'Cdc25i 10uM + 4/6i',5,2:6,1:4;
%     %'Wee1i 100nM',6,2:6,1:4;
%     %'Wee1i 100nM + 4/6i',7,2:6,1:4;
%     };

%%% 20131213 R-point CC\20140506 3C G1S TranxTrans
conditions={
    'sample',1:2,1,1;
    };

%%% GENERAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='H:\Documents\Projects\';
%experimentpath='20131213 R-point CC\20131217 3Doldrabbit XPcdkpanel\';
%experimentpath='20131126 R-point SR\20140317 3C CDK4i Wee1i\';
%experimentpath='20131126 R-point SR\20140323 Wee1i Cdc25i Cdk2i\';
%experimentpath='20131126 R-point SR\20140324 3C CDK4i CDK2i\';
%experimentpath='Heewon\';
%experimentpath='20131209 Pipdeg Optimization\20140402 PFACS\';
%experimentpath='20131213 R-point CC\20140318 CyclingCDK4iWee1iCDC25i\';
%experimentpath='20131213 R-point CC\20140406 siRNA CDK4i\';
%experimentpath='20131126 R-point SR\20140404 siCycE siCycD CDK4i CDK2i\';
%experimentpath='20131213 R-point CC\20140407 2C CDK2 Hysteresis\';
%experimentpath='20131213 R-point CC\20140413 2C DHB Hysteresis\';
%experimentpath='20131126 R-point SR\20140415 CDK46hysteresis-under-CDK2i\';
%experimentpath='20131209 Pipdeg Optimization\20140405 PFACS\';
%experimentpath='Heewon\BreakBeforeMitosis\';
%experimentpath='20131213 R-point CC\20140423 3C G1S 46i\';
experimentpath='20131213 R-point CC\20140505 3C G1S TranxTrans\';

datadir=([projectpath,experimentpath,'Data\']);
%datadir=([projectpath,experimentpath,'Data_ksmode\']);
%datadir=([projectpath,experimentpath,'Data_bgjustmean\']);
%datadir=([projectpath,experimentpath,'Data_ConstantBG_0th_min\']);
%datadir=([projectpath,experimentpath,'Data_min\']);
%datadir=([projectpath,experimentpath,'Data_min_ringexcludebg\']);
%datadir=([projectpath,experimentpath,'Data_min_ringexcludebg_rad4\']);

%%% CHOOSE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CellCount(conditions,rawpath,nucr,1,[1 226]);
%IMT(conditions,datadir);
%DAPI_vs_EdU(conditions,datadir,1);
%DAPI_vs_sensor(conditions,datadir);
%DAPI_vs_DHB(conditions,datadir);
%DHB_vs_sensor(conditions,datadir);
%sensor_vs_EdU_test(conditions,datadir);
%Align_Drugspike(conditions,datadir);
%Timelapse_IF(conditions,datadir);
%PhaseLengths_genealogy(conditions,datadir);
%Cdk2(conditions,datadir);
%Cdk2_pipdeg(conditions,datadir);
%Timelapse_IF_Cdk2_pipdeg(conditions,datadir);
%Geminin_Comparison(conditions,datadir);
%Geminin_p21(conditions,datadir);

%%% validated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ShadingEval(conditions,datadir);
%Stain_TrendComparison(conditions,datadir);
%Stain_HistComparison(conditions,datadir);
%Stain_3Dmap(conditions,datadir);
%Geminin(conditions,datadir);
%pipdeg_CDK2(conditions,datadir);
%Timelapse_IF_Cdk2(conditions,datadir);
%Visualize_CellCycleGating(conditions,datadir);
Timelapse_General_G1S(conditions,datadir);
%Timelapse_General_CellCycle(conditions,datadir);
%Timelapse_pRb_vs_R(conditions,datadir);