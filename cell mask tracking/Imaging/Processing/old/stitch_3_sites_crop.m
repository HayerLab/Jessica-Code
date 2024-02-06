%% Stitch together 3 vertically tiled sites 
clear;clc;
close all;

root='N:\161103_IXM\raw';
datadir=[root,filesep,'stitched'];
if ~exist(datadir)
    mkdir(datadir)
end

for row=1:8
    for col=1:12
        position=[num2str(row),'_',num2str(col),'_'];
        for frame=[1 50:70];
            temp1=imread([root,filesep,position,'1_DAPI_',num2str(frame),'.tif']);
            temp2=imread([root,filesep,position,'2_DAPI_',num2str(frame),'.tif']);
            temp3=imread([root,filesep,position,'3_DAPI_',num2str(frame),'.tif']);

            numrows=size(temp1,1);
            stitched=zeros(3*numrows,size(temp1,2));

            stitched=temp1;
            stitched((numrows+1):2*numrows,:)=temp2;
            stitched((2*numrows+1):3*numrows,:)=temp3;
            %crop
            stitched=stitched(1621:4860,:);
            
            
            imwrite(stitched, [datadir,filesep,position,'1_DAPI_',num2str(frame),'.tif']);
            disp(position);
        end
    end
end
%imagesc(stitched);
disp('done!');

