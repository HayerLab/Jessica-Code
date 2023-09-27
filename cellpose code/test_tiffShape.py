# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 11:01:44 2023

@author: jessicaz
"""

import tifffile
import os.path
filename = os.path.join( "segmentation", "3_7_1", "Composites","aligned","composite_aligned_1.tif" )
data = tifffile.imread(filename)
print(data.shape)

example_path = os.path.join("Analyzing_fluorescence_images_data",
                            "Analyzing_fluorescence_images_data","Cell_composite.tif")
data = tifffile.imread(example_path)
print(data.shape)