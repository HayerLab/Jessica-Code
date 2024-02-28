# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 19:35:11 2024

@author: zhuje
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 17:31:58 2024

@author: zhuje
"""

import pandas as pd
import scipy as scp
import numpy as np
import seaborn
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error,mean_squared_error,accuracy_score,r2_score
import sys
sys.path.append('Z:\\Jessica\\tracking_code\\factorial_hmm')
import factorial_hmm
from filterpy.kalman import KalmanFilter
import matplotlib.pyplot as plt

#load MATLAB file
data = scp.io.loadmat('Z:\\Jessica\\tracking_code\\data\\varied_model\\3_8_1_mask_tracking_varied_model')
filtered = scp.io.loadmat('Z:\\Jessica\\tracking_code\\label_matrix\\single_nuclei\\data\\3_8_1\\labels\\90_filtered_labels_consistent_3_8_1.mat')
cell_filter = filtered.get('common_val')
filter_list = np.ndarray.tolist(cell_filter.reshape(1,-1)[0])
mask_data = data.get('mask_track')

#relevant cell shape and motility parameter values
v = mask_data[:,:,5]
v = np.nan_to_num(v)
area = mask_data[:,:,12]
perim = mask_data[:,:,17]
ecc = mask_data[:,:,16]
circ = mask_data[:,:,15]
area = np.nan_to_num(area)
perim = np.nan_to_num(perim)
ecc = np.nan_to_num(ecc)
circ = np.nan_to_num(circ)
mean_ecc = np.mean(ecc,1)
mean_circ = np.mean(circ,1)

# 1 pixel ~ 0.325um
v = v * 0.325
area = area * 0.325 **2
perim = area * 0.325

#cell shape index defined as P/sqrt(A)
cell_shape = np.mean(perim,1)/np.sqrt(np.mean(area,1))
nan_idx = np.argwhere(~np.isnan(cell_shape))
mean_v = np.mean(v,1)
mean_v = mean_v[nan_idx]

cell_shape = cell_shape[nan_idx]
mean_ecc = mean_ecc[nan_idx]
mean_circ = mean_circ[nan_idx]

#linear regression model with average cell shape parameters as features
features = np.hstack((cell_shape[filter_list],mean_ecc[filter_list]))
X_train, X_test, y_train, y_test = train_test_split(features,mean_v[filter_list],test_size = 0.20, train_size = 0.80, random_state = 32)
model = LinearRegression().fit(X_train,y_train)
predictions = model.predict(X_test)
metric = model.score(X_test,y_test)

all_shape = perim/np.sqrt(area)
all_shape = np.nan_to_num(all_shape)
df_shape = pd.DataFrame(all_shape[filter_list,0:179])
mean_shape_seq = df_shape.mean(axis=0)

#distributions
fig, ax = plt.subplots(figsize=(15, 6))
seaborn.distplot(mean_shape_seq, kde = True, color = 'blue', ax = ax)
ax.set_title("Distribution of Cell Shape")
ax.set_xlabel("Cell Shape Index")
df_v = pd.DataFrame(v[filter_list,0:179])
correlation = df_shape.corrwith(df_v,axis=0)
mean_v_seq = df_v.mean(axis=0)
fig_v, ax_v = plt.subplots(figsize=(15, 6))
seaborn.distplot(mean_v_seq, kde = True, color = 'blue', ax = ax_v)
ax_v.set_title("Distribution of Velocity")
ax_v.set_xlabel('Velocity')

meanVDF = pd.DataFrame(mean_v)
meanShapeDf = pd.DataFrame(np.mean(circ,1))
mean_corr = meanVDF.corrwith(meanShapeDf,axis=0)
df_ecc = pd.DataFrame(ecc[filter_list,0:179])
df_circ = pd.DataFrame(circ[filter_list,0:179])
df_pred = pd.DataFrame()
for row in range(len(filter_list)):
    for col in range(0,179):
        df_pred.loc[row,col] = 0


for cell in range(len(filter_list)):
    for f in range(0,179):
        #linear predictor based on last 2 datapts
        if f > 1:
            prev = df_v.loc[cell,f-1]
            next_prev = df_v.loc[cell,f-2]
            rise =  prev - next_prev
            shape_prev = np.asarray([df_shape.loc[cell,f-1],df_ecc.loc[cell,f-1],df_circ.loc[cell,f-1]])
            shape_next = np.asarray([df_shape.loc[cell].values[f],df_ecc.loc[cell].values[f],df_circ.loc[cell].values[f]])
            X_sample = shape_prev.reshape(1,-1)
            prev_pred = model.predict(X_sample)
            next_pred = model.predict(shape_next.reshape(1,-1))
            pred_diff = next_pred - prev_pred
            change_diff = pred_diff - rise
            df_pred.loc[cell,f] = next_pred + abs(change_diff)
                 
variance = np.var(v[filter_list,0:179],axis=0)
actual = df_v.loc[:,2:179]
predictions = df_pred.loc[:,2:179]
diff = np.subtract(actual,predictions)              
mse = mean_squared_error(actual,predictions)
mae = mean_absolute_error(actual, predictions)
r2 = r2_score(actual, predictions)
mean_var = np.mean(variance)
plt.plot(list(range(2,179)),actual.loc[62,:],label='Observed')
plt.plot(list(range(2,179)),predictions.loc[62,:], label='Predicted')
plt.xlabel('Frame')
plt.ylabel('Velocity (um/min)')
plt.legend() 


plt.scatter(df_shape.mean(axis=1),df_v.mean(axis=1),marker='.', s=15)
plt.xlabel('Cell Shape Index (um)')
plt.ylabel('Single Cell Velocity (um/min)')
plt.show()


sampled = np.random.choice(all_shape.shape[1], 20, replace=False)
X = all_shape[:, sampled]
y = v[:, sampled]
test_X = np.delete(all_shape,sampled, axis=1)
test_y = np.delete(v, sampled,axis=1)
single_X = X[95,:]
single_y = y[95,:]





