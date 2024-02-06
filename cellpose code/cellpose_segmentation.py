# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 13:43:05 2023

@author: jessicaz
"""
import numpy as np
import time, os, sys
from urllib.parse import urlparse
import skimage.io
import matplotlib.pyplot as plt
import matplotlib as mpl
import cellpose
from urllib.parse import urlparse
from cellpose import models, core
from tifffile import imread, imwrite
import torch
mpl.rcParams['figure.dpi'] = 300
use_GPU = core.use_gpu()
from cellpose.io import logger_setup
import os
logger_setup();
# setting device on GPU if available, else CPU
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print('Using device:', device)
print()

#Additional Info when using cuda
if device.type == 'cuda':
    print(torch.cuda.get_device_name(0))
    print('Memory Usage:')
    print('Allocated:', round(torch.cuda.memory_allocated(0)/1024**3,1), 'GB')
    print('Cached:   ', round(torch.cuda.memory_reserved(0)/1024**3,1), 'GB')

drive = "Z:\\"
user = "Jessica\\"
base_folder = "segment\\"
shots =["3_8_1\\", "3_8_2\\"]
model_filepath = "model with varied sets\\training_set\\models\\HUVEC_monolayer_varied_set"
image_folder = "Composites\\aligned\\tif\\" #note: 2 backward slash for windows folder, 1 for mac
#model_path = "E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig_tiff\\training\\models\\Improved_HUVEC_monolayer" #@param {type:"string"}
model_path = os.path.join(drive, user, model_filepath)
#@markdown ###Path to images:

for shot in shots:
#directory = "E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned" #@param {type:"string"}
    directory = os.path.join(drive, user, base_folder, shot, image_folder)
    #@markdown ###Channel Parameters:
    
    Channel_to_use_for_segmentation = "Green" #@param ["Grayscale", "Blue", "Green", "Red"]
    
    # @markdown If you have a secondary channel that can be used, for instance nuclei, choose it here:
    
    Second_segmentation_channel= "Red" #@param ["None", "Blue", "Green", "Red"]
    
    
    # Here we match the channel to number
    if Channel_to_use_for_segmentation == "Grayscale":
      chan = 0
    elif Channel_to_use_for_segmentation == "Blue":
      chan = 3
    elif Channel_to_use_for_segmentation == "Green":
      chan = 2
    elif Channel_to_use_for_segmentation == "Red":
      chan = 1
    
    
    if Second_segmentation_channel == "Blue":
      chan2 = 3
    elif Second_segmentation_channel == "Green":
      chan2 = 2
    elif Second_segmentation_channel == "Red":
      chan2 = 1
    elif Second_segmentation_channel == "None":
      chan2 = 0
    
    #@markdown ### Segmentation parameters:
    
    #@markdown diameter of cells (set to zero to use diameter from training set):
    diameter =  105#@param {type:"number"}
    #@markdown threshold on flow error to accept a mask (set higher to get more cells, e.g. in range from (0.1, 3.0), OR set to 0.0 to turn off so no cells discarded):
    flow_threshold = 0.4 #@param {type:"slider", min:0.0, max:3.0, step:0.1}
    #@markdown threshold on cellprob output to seed cell masks (set lower to include more pixels or higher to include fewer, e.g. in range from (-6, 6)):
    cellprob_threshold=0 #@param {type:"slider", min:-6, max:6, step:1}
    
    # gets image files in directory (ignoring image files ending in _masks)
    files = cellpose.io.get_image_files(directory, ('_seg'))
    print(files)
    images = [cellpose.io.imread(f) for f in files]
    #images = cellpose.io.imread("E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\composite_aligned_58.tif")
    # declare model
    print("run")
    model = models.CellposeModel(gpu=True, device=torch.device('cuda'),
                                 pretrained_model=model_path)
    
    print("eval")
    # run model on test images
    masks, flows, styles = model.eval(images, 
                                      channels=[chan, chan2],
                                      diameter=105,
                                      flow_threshold=flow_threshold,
                                      cellprob_threshold=cellprob_threshold
                                      )
        
    # cellpose.io.save_masks(images, 
    #               masks, 
    #               flows, 
    #               files, 
    #               channels=[chan, chan2],
    #               png=True, # save masks as PNGs and save example image
    #               tif=True, # save masks as TIFFs
    #               save_txt=True, # save txt outlines for ImageJ
    #               save_flows=False, # save flows as TIFFs
    #               save_outlines=False, # save outlines as TIFFs 
    #               )
        
    
    cellpose.io.masks_flows_to_seg(images, 
                          masks, 
                          flows, 
                          diameter*np.ones(len(masks)), 
                          files, 
                          [chan, chan2])
    
    # from cellpose import plot
    # channels = [2,1]
    # nimg = len(files)
    # for idx in range(nimg):
    #     maski = masks[idx]
    #     flowi = flows[idx][0]
    #     fig = plt.figure(figsize=(12,5))
    #     plot.show_segmentation(fig, images[idx], maski, flowi, channels=[2,1])
    #     print(idx)
    #     plt.tight_layout()
    #     plt.show()
    #   # plt.savefig("C:\\Users\\Jessica Zhu\\Downloads\\segmentation\\3_7_1\\Composites\\aligned\\export_fig tiff\\results\\test_set.png")
        
    i = 0
    im = "images\\"
    csv = "csv\\"
    newpath_image = os.path.join(drive, user, base_folder, shot, image_folder, im)
    newpath_csv = os.path.join(drive, user, base_folder, shot, image_folder, csv)
    if not os.path.exists(newpath_image):
        os.makedirs(newpath_image)
    if not os.path.exists(newpath_csv):
        os.makedirs(newpath_csv)
    for mask in masks: 
        
            #imwrite(directory+"\\images\\"+str(i+1)+".tif", mask) 
        imwrite(os.path.join(newpath_image, str(i+1)+".tif"), mask)
        #np.savetxt(directory+"\\csv\\segmented_"+str(i+1)+".csv", mask, delimiter=",")
        np.savetxt(os.path.join(newpath_image, "segmented"+str(i+1)+".csv"), mask, delimiter=",")
        i+=1
        
    print("save masks")
    
    
