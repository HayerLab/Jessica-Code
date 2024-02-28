# Whole-cell Segmentation and Tracking for HUVEC Monolayers
This repository contains the codebase for whole-cell segmentation and tracking for HUVEC monolayers. Segmentation is done using Cellpose and cell tracking is done using a series of Matlab scripts. 

Contributors: Jessica, Lydia

## Table of contents
1. [Preprocessing](#preprocessing)
2. [Cell Segmentation](#cell-segmentation)
    1. [Pre-requisite](#pre-requisite)
    2. [General Workflow](#general-workflow)
3. [Cell Tracking](#cell-tracking)
    1. [Matching Cellpose Masks to Nuclear Tracking Output](#matching-cellpose-masks-to-nuclear-tracking-output)
        1. [Setup](#setup)
        2. [Running the script](#running-the-script)
    2. [Postprocessing](#postprocessing)
        1. [Cell number labeling](#cell-number-labeling)
        2. [Remove border cells](#remove-border-cells)
        3. [Intensity filtering](#intensity-filtering)
        4. [Cell selection based on track length](#cell-selection-based-on-track-length)


## Preprocessing 
The raw data consists of 16-bit greyscale fluorescence microscopy images of the nuclear and E-cadherin stains. Preprocessing is done to generate RGB images of the monolayer from the raw images with DAPI as the red channel and YFP as the green channel. The Preprocessing folder contains code for generating the composite images, as well as the required dependencies. 

First, open the ``` GenerateReferenceBackgroundImages.m ``` file. Change line 4 and 5 according to your directory names: 
```
rawpath='<folder containing raw background images>'; 
background='<name of background folder>';
```
Now run the file. 

Next, open ``` preprocessing/AlignmentParameters.m ```. Change the shot variable to the shot number you are running on and change the root, rawdir and bgdir variables according to your folder names. 

Then, open ``` preprocessing/GenerateCompositeImages.m ```. Change the root, rawdir, bgdir and merged_channel_dir variables to the corresponding folder names on your computer.

Example: 
```
root = 'Z:\Jessica';
rawdir=[root,filesep,'raw',filesep,'160519_CDH5_mCitrine_Hoechst_1min_interval_20x',filesep, '3_7_1'];
bgdir=[root,filesep,'background'];
merged_channel_dir = [root, filesep,'raw', filesep, 'composites', filesep, '3_7_1'];
```

 The ```load ``` command in line 21 loads the stored Matlab matrix output from ```AlignmentParameters.m```, change the filepath accordingly. The position variable should also be changed to the shot number you are currently analyzing. 

## Cell Segmentation

### Pre-requisite
You must have the Cellpose library installed prior to running any cell segmentation code. For details on how to install Cellpose, please see [Cellpose GitHub Repository](https://github.com/MouseLand/cellpose).


### General Workflow
Cell segmentation is done using a custom trained Cellpose model. A python script is setup to automate the segmentation process. Python scripts for training a cellpose model, as well as analyzing its performance are also stored in the ```cellpose code```  folder. 

There are two HUVEC-monolayer-trained cellpose models available, both under the ```trained cellpose models``` folder. One is trained and tested on a single dataset consisting of images from a single movie. The other is trained with a dataset that consists of images from 3 different movies (15 images each). The two models have similar performances on segmenting HUVEC monolayers, so it doesn't matter which one you choose to use. 

To perform cell segmentation, open the ```cellpose_segmentation.py``` file. The following path variables need to be changed accordingly: 
```
# note: all folder names must end with a file separator (\\ for windows)
drive = "<base drive folder>"
user = "<username folder>"
base_folder = "<base folder name containing folder with composite images>"
shots =[<list of the shot number as strings>]
model_filepath = "<filepath for cellpose model>"
image_folder = "<filepath to image folder>"
```
Now you can run the file. The segmentation results are stored at np files in the image folder, but also stored as .tif images and .csv files under separate images and csv folders within the original image folder. 
The .csv files contain the label matrices of the segmented cell masks and we will use these files in later cell tracking. 

## Cell Tracking 

### Matching Cellpose Masks to Nuclear Tracking Output

The top level script for cell tracking is the ``` mask_matching.m ``` file. The cell mask tracking consists of three parts: nuclear tracking, motility parameters computations, cell mask matching. 

#### Setup
Change the ``` shots_row```, ``` shots_col```, and ``` shots_sites``` variables on line 25 accordingly. The contents of the cell array depends on how many and what datasets are you running the tracking code on

Line 3 to 9 sets up the folder for the segmented masks, raw images, data storage etc. Please change these accordingly.

Depending on the dataset that you are running the tracking algorithm on, you will have to adjust the input parameters to the nuclear tracking code. The Timelapse20x() function takes as input the respective row,col,site of the microscopy images,projectpath,imagepath,experimentpath, SF (starting frame), and EF (ending frame). See below for an example usage for a movie with 25 frames:
```
nuc_tracking = Timelapse20x(shot_num(1), shot_num(2),shot_num(3),projectpath,imagepath,experimentpath, 1,25);

```
#### Running the script
The first two parts are done with previous cell tracking scripts that are contained within the Imaging and MotilityParameters folders. When the mask_matching script is run, it automatically run ```add_path.m``` which sets up the correct paths. 

The output of running the nuclear tracking and motility parameters code is an N x M x 11 3-D matrix where N is the number of cells and M is the number of frames. This matrix would be saved as a MAT-file object. 

The function ```mask_matching()``` is called to perform the mask matching. The function takes in the matrix output from the previous two steps and outputs an N x M x 18 3-D matrix. The matrix is saved as a MAT-file object. 

Each slice of the matrix contains the values for different parameters of each cell. The parmeters are as follow:

1. x-coordinates of nulear centroids
2. y-coordinates of nuclear centroids
3. Area of nucleus
4. Mass of nucleus (density x area)
5. Angle of movement 
6. Single cell velocity
7. Coordination with neighboring cells in the front, located within a 60° sector ahead
8. Coordination with lateral neighbors, located within lateral 120° sectors
9. Coordination with neighboring cells in the back, located within a 60° sector behind
10. Coordination with all neighboring cells
11. Persistence. col 1: net distance. Col2: total distance. Col3: net/total distance.
12. Matched cellpose mask labels
13. Area of cell masks
14. x-coordinates of cell mask centroids
15. y-coordinates of cell mask centroids
16. Circularity of cell masks
17. Eccentricity of cell masks
18. Perimeter of cell masks

### Postprocessing
The output of the cell mask tracking is postprocessed for better downstream analysis. Postprocessing is done in three steps: create new label matrices with cell number for single nucleu cells only, remove cells at the border, filter out cells with higher than average intensity. 

#### Cell number labeling
Label matrices (same size as original label matrix output from the cell segmentation) are generated for each frames with each label representing the assigned cell number. As part of the relabeling, cells with multiple nuclei are excluded from the label matrix for the purpose of downstream analyses. However, the code can be modified to include such cells. 

The top level script is ```cell_number_labels.m```. The script calls the ```cell_number``` function (see cell_number.m) which generates the new label matrix. For the purpose of genrating and saving label matrices of cell numbers, the function returns a M x N size 2-D matrix where M, N are the size of the original label matrix. 

The current script excludes all non-single nucleus cells when labeling the matrix. If all cells (single nucleus or not) should be included then replace line 21 - 28: 
```
if isnan(cell_label(label_idx)) & cell_label(label_idx)~= 0
    ID = cellID(i);
    cell_label(label_idx) = ID;
elseif ~isnan(cell_label(label_idx))
    cell_label(label_idx) = 0;
elseif cell_label(label_idx) == 0
    continue;
end
```
with 
```
ID = cellID(i);
cell_label(label_idx) = ID;
```

#### Remove border cells
Cells at or near the border of the image will be removed since these cells are generally only partially in the field of view (part of its cell is cut off). The main script for this process is ```remove_border.m```. 

The following lines of code must be changed depending on the dataset used: 
- Line 10 - 14: change the file paths and rows/cols/sites variables accordingly
- Line 55: change the range to the number of frames for the dataset you are analyzing

You can also change the border threshold input of the good_labels() function to adjust the subset of cells that are at least x pixels away from the image border you want to keep. 

The code saves all the new cell number matrices for each frame and also saves a struct array containing the filtered cells. 

#### Intensity filtering
Although the trained cell segmentation model generally correctly segments cells, it does have difficulty segmenting individual cells in a "crowded" region of the image (possibly from multiple cells whose edges overlap with each other). Such cells are removed from the label matrix by filtering out any cells that has a high mean intensity. The threshold used for the filtering is the averagee of all computed cell mean intensity plus one standard deviation. 

Similar to the previous step, path variable, as well as row/col/site variables should be changed accordingly. 

Line 39 and 67 should also be changed to match the range with the number of frames in the dataset being used.

The code saves all the new cell number matrices for each frame and also saves a struct array containing the filtered cells (i.e. cells whose mena intensity did not exceed the threshold). 

#### Cell selection based on track length
This step filters out any cells whose track length is less than some set threshold. 

To adjust the threshold, change the threshold variable declared in Line 14 in the consistent_tracks.m file. 

As before, path variables etc. should also be changed accordingly. Change Line 48 such as the number of for loop iterations matches the number of frames in your dataset. 

