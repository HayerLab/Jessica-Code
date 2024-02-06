function FormatFiles(row,col,site)
%row=3; col=5; site=1;
stain=0; %1:immunostain 0:timelapse
separatedirectories=0;
time2=tic;

rawdir='K:\160831-AH\'; % Place where the new files will be

[rowstring,colstring,sitestring]=wellnum2strRCS_2(row,col,site); %convert integers to A_01_1 format
shot=[num2str(row),'_', num2str(col), '_', num2str(site)];

camdir='K:\160831-AH\2016-08-31\61'; % files where the old files are

w1root='DAPI';
% w2root='YFP';
% w3root='RFP';
% w4root='Cy5';
%w5root='EdU';

drugtime=0;%33; %0,39,103
times=1:25;

shotdir=[rawdir,'\',shot];
if ~separatedirectories
    shotdir=rawdir;
    w1root=[shot,'_',w1root];
%     w2root=[shot,'_',w2root];
%     w3root=[shot,'_',w3root];
%     w4root=[shot,'_',w4root];
    %w5root=[shot,'_',w5root];
elseif ~exist(shotdir,'dir')
    mkdir(shotdir);
end
shotdir=[shotdir,'\'];
for time=times
    if stain
        timedir=camdir; %Immunostain
    else
        timedir=[camdir,'\TimePoint_',num2str(time)]; %Timelapse
    end
    filenames=getFilenames(timedir);
    filenames=filenames(~boolRegExp(filenames,'thumb'));
    
    %%% in case of only one channel
    %%%w1old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring,'_s',sitestring]))); %for when only one channel is captured
    
    %%% in case of only one site
    %w1old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring,'_w1'])));
    %w2old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring,'_w2'])));
    %w3old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring,'_w3'])));
    %w4old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring,'_w4'])));
    
    %%% multiple sites and channels
 w1old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring])));
    %w1old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring,'_s',sitestring,'_w1'])));
    %w2old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring])));
    %w3old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring])));
%     w4old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring,'_s',sitestring,'_w4'])));
    %w5old=char(filenames(boolRegExp(filenames,['_',rowstring,colstring,'_s',sitestring,'_w5'])));
    
    if stain
        newstring='stain';
        %newstring='shading';
        %newstring='serum';
    else
        newstring=num2str(time+drugtime);
    end
    w1new=[w1root,'_',newstring,'.tif'];
    
%     w2new=[w2root,'_',newstring,'.tif'];
%     w3new=[w3root,'_',newstring,'.tif'];
%     w4new=[w4root,'_',newstring,'.tif'];
    %w5new=[w5root,'_',newstring,'.tif'];
    %system(['copy "',timedir,'\',w1old,'" "',shotdir,w1new,'"']);
    %system(['copy "',timedir,'\',w2old,'" "',shotdir,w2new,'"']);
    %system(['copy "',timedir,'\',w3old,'" "',shotdir,w3new,'"']);
    %system(['copy "',timedir,'\',w4old,'" "',shotdir,w4new,'"']);
    %system(['copy "',timedir,'\',w5old,'" "',shotdir,w5new,'"']);
    system(['move "',timedir,'\',w1old,'" "',shotdir,w1new,'"']);
    %system(['move "',timedir,'\',w2old,'" "',shotdir,w2new,'"']);
    %system(['move "',timedir,'\',w3old,'" "',shotdir,w3new,'"']);
%     system(['move "',timedir,'\',w4old,'" "',shotdir,w4new,'"']);
    %system(['move "',timedir,'\',w5old,'" "',shotdir,w5new,'"']);
end
        
%%% reference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%system(['mkdir ',path,'Raw']);
%system(['move ',path,'testb.tif ',path,'Raw']);
%system(['ren ',path,'testf.tif testg.tif']);
%system(['copy ',path,'testa.tif ',path,'Raw\testacopied.tif']);
%system(['del ',fldir,'* /Q']);
toc(time2)