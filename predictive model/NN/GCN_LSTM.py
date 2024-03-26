#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 10:26:43 2024

@author: jessicazhu
"""

import os
import sys
import urllib.request

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy as scp
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import Sequential, Model
from tensorflow.keras.layers import LSTM, Dense, Dropout, Input
import pandas as pd

import stellargraph as sg
from stellargraph.layer import GCN_LSTM
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))


def train_test_split(data, train_portion):
    time_len = data.shape[1]
    train_size = int(time_len * train_portion)
    train_data = np.array(data.iloc[:, :train_size])
    test_data = np.array(data.iloc[:, train_size:])
    return train_data, test_data

data = scp.io.loadmat('Z:\\Jessica\\tracking_code\\data\\updated\\3_8_1_mask_tracking_varied_v2')
filtered = scp.io.loadmat('Z:\\Jessica\\tracking_code\\label_matrix\\single_nuclei\\data\\3_8_1\\labels\\90_filtered_labels_consistent_3_8_1.mat')
cell_filter = filtered.get('common_val')
filter_list = np.ndarray.tolist(cell_filter.reshape(1,-1)[0])
mask_data = data.get('mask_track')
v = mask_data[:,:,5]

area = mask_data[:,:,14]
perim = mask_data[:,:,18]
ecc = mask_data[:,:,17]
circ = mask_data[:,:,16]
coord_front = mask_data[:,:,6]
coord_lateral = mask_data[:,:,7]
coord_back = mask_data[:,:,8]
neighbours = mask_data[:,:,11]
coord_nb = mask_data[:,:,9]
v = v * 0.325
area = area * 0.325 **2
perim = perim * 0.325
cell_shape = perim/np.sqrt(area)
shape_filtered = cell_shape[filter_list]
NbCoordFiltered = coord_nb[filter_list]
df_v = pd.DataFrame(v)
df_shape = pd.DataFrame(cell_shape)
adj_data = scp.io.loadmat('Z:\\Jessica\\tracking_code\\data\\adjacency\\3_8_1\\adjacency_flattened_3_8_1')
adj_mat = adj_data.get('adj_mat')

train_rate = 0.8
train_data, test_data = train_test_split(df_shape, train_rate)
print("Train data: ", train_data.shape)
print("Test data: ", test_data.shape)

def scale_data(train_data, test_data):
    max_val = train_data.max()
    min_val = train_data.min()
    train_scaled = (train_data - min_val) / (max_val - min_val)
    test_scaled = (test_data - min_val) / (max_val - min_val)
    return train_scaled, test_scaled
train_scaled, test_scaled = scale_data(train_data, test_data)
trainX = train_data
testX = test_data



trainY, testY =  train_test_split(df_v, train_rate)
trainY,testY = scale_data(trainY,testY)



print(trainX.shape)
print(trainY.shape)
print(testX.shape)
print(testY.shape)
seq_len = 10
pre_len = 12
gcn_lstm = GCN_LSTM(
    seq_len=seq_len,
    adj=adj_mat,
    gc_layer_sizes=[16, 10],
    gc_activations=["relu", "relu"],
    lstm_layer_sizes=[200, 200],
    lstm_activations=["tanh", "tanh"],
)

x_input, x_output = gcn_lstm.in_out_tensors()
model = Model(inputs=x_input, outputs=x_output)

model.compile(optimizer="adam", loss="mae", metrics=["mse"])
history = model.fit(
    trainX,
    trainY,
    epochs=100,
    batch_size=60,
    shuffle=True,
    verbose=1,
    validation_data=[testX, testY],
)

model.summary()




