
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:20:16 2023
modified from https://github.com/MouseLand/cellpose/issues/496
move masks to one directory using zsh command cp **/*.npy ./all_masks
"""
import pickle
import pandas as pd
import scipy
import numpy as np
import cellpose
import natsort
from os import listdir
import matplotlib.pyplot as plt
from datetime import datetime
import gc

print(datetime.now())

# Function for calculating centerpoints (again)
def calculate_centers(masks):
    slices = scipy.ndimage.find_objects(masks)
    centroids = list()
    centroids_x = list()
    centroids_y = list()
    for i, elem in enumerate(slices):
        if elem is not None:
            sr, sc = elem
            mask = masks[sr, sc]
            center = scipy.ndimage.center_of_mass(mask)
            centroids.append(center)
    return centroids
#directory = ["E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig_tiff\\all_masks\\sub\\"]
directory = ["E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig_tiff\\all_masks\\","E:\\Jessica\\segmentation\\3_8_2\\Composites\\aligned\\"]
shots = ["3_7_1","3_8_1", "3_8_2"]
# "E:\\Jessica\\segmentation\\3_8_1\\Composites\\aligned\\", 

# function for loading information and calculate center positions
def get_centroid(directory):
    gc.collect() #garbage collection to free up memory
    data =dict()
    masks = dict()
    centers = dict()
    image = dict()
   
    #load data
    for shot in directory:
        print(shot)
        k = 0
        shot_data = list()
        sorted_dir_content = natsort.natsorted(listdir(shot))
        for content in sorted_dir_content:
            if content.split('.')[-1] == "npy":
                print(k)
                #load and append data for each shot
                shot_data.append(np.load(shot+content, allow_pickle=True).item())
                gc.collect()
                k+=1
        data[shot] = shot_data
    frames_mask = dict()
    frames_img = dict()
    value = data.values()
   
    #get mask and images for each shot
    for key,value in data.items():
            masks_timeseries = list()
            image_series = list()
            for timepoint in value:
                masks_timeseries.append(timepoint['masks'])
                if key != 'E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig_tiff\\all_masks\\':
                        image_series.append(timepoint['img'])
                masks[key] = masks_timeseries
                image[key] = image_series

    
    #calculate center for each frame for each shot 
    center_df_shots = dict()

    for key,value in masks.items(): #per shot
        center_df_frames = dict()
        center_timepoint_x = list()
        center_timepoint_y = list()
        columns = list(range(len(value)))
        df =  pd.DataFrame(columns=columns)
        df_y = pd.DataFrame(columns=columns)
        for index, elem in enumerate(value): #per frame
            center = calculate_centers(elem)
            df= pd.DataFrame(center)
            center_df_frames[index] = df
        center_df_shots[key] = center_df_frames
    return data,masks,image,center_df_shots

data,masks,image,centroid = get_centroid(directory)
print(datetime.now())
index = 0
for k,v in centroid.items(): 
    for i,frame in v.items():
        frame.to_csv("E:\\Jessica\\segmentation\\centroid\\Centroid"+shots[index]+"_"+str(i)+".csv")
    index+=1
