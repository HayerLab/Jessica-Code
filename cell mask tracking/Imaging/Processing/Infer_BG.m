function Average_ShadingImages
%%% file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagepath='F:\';
%imagepath='E:\';
shadingpath='H:\Images\ShadingImages\20140425 DCYTC 10x\';
%shadingpath='H:\Images\ShadingImages\20140522 20xBin2\';
experimentpath='Michael-20140608\Start\2014-06-08\';
separatedirectories=1;
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names={
    'CFP';
    'YFP';
	%'tRb';
	%'pRb';
	%'Cy5';
};
rowmat=[3:6]; %5 6
colmat=[2:11]; %5 6
sitemat=1:4; %1
framemat=[1];
nucr=12;
%%% average images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows=numel(rowmat); numcols=numel(colmat); numsites=numel(sitemat); numframes=numel(framemat);
load([shadingpath,'BG.mat'],'shadingcorrection'); bgcmos=shadingcorrection;
for i=1:length(names)
    name=char(names{i});
    for s=1:numsites
        rawarrays=[];
        site=sitemat(s);
        for c=1:numcols
            col=colmat(c);
            for r=1:numrows
                row=rowmat(r);
                for f=1:numframes
                    frame=framemat(f);
                    wellcolname=wellcol(row,col);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    if separatedirectories
                        %rawdir=[imagepath,experimentpath,'Raw\'];
                        rawdir=[imagepath,experimentpath,'Raw\',wellcolname,'\site_',num2str(site),'\'];
                        %raw=single(imread([rawdir,names{i},'_',num2str(frame),'.tif']));
                        raw=single(imread([rawdir,shot,'_',names{i},'_',num2str(frame),'.tif']));
                    else
                        rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
                        raw=single(imread([rawdir,names{i},'_stain.tif']));
                    end
                    raw=raw-bgcmos;
                    rawarrays=cat(3,rawarrays,raw);
                end
            end
        end
        rawcalc=min(rawarrays,[],3);
        inferredbg=shadingcorrection_1([],rawcalc,3*nucr); %can be 1*nucr with enough sites
        save([imagepath,experimentpath,'Raw\IBG\',names{i},'_',num2str(site),'.mat'],'inferredbg');
        %save([imagepath,experimentpath,'Raw\IBG_',names{i},'.mat'],'inferredbg');
    end
end