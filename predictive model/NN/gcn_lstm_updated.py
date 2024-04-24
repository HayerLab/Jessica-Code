# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:01:46 2024

@author: zhuje
"""


"""
## Setup
"""

import os

os.environ["KERAS_BACKEND"] = "tensorflow"
import random
import pandas as pd
import numpy as np
import seaborn as sns
import typing
import matplotlib.pyplot as plt
from random import sample
import tensorflow as tf
print(tf.__version__)
import keras
from keras import layers
from keras import ops
import scipy as scp
from keras import backend as K
from sklearn.metrics import classification_report


"""
### Loading data
"""
frame = 180
data_short = scp.io.loadmat('/mnt/d/Jessica/2_4_2_mask_tracking.mat')
masks = data_short.get('mask_track')
data = scp.io.loadmat('/mnt/d/Jessica/3_8_1_mask_tracking_varied_model.mat')
filtered = scp.io.loadmat('/mnt/d/Jessica/180_filtered_labels_consistent_3_8_1.mat')
cell_filter = filtered.get('common_val')
filter_list = np.ndarray.tolist(cell_filter.reshape(1,-1)[0])
mask_data = data.get('mask_track')
data2 = scp.io.loadmat('/mnt/d/Jessica/3_8_2_mask_tracking_varied_model.mat')
filtered2 = scp.io.loadmat('/mnt/d/Jessica/180_filtered_labels_consistent_3_8_2.mat')
cell_filter2 = filtered2.get('common_val')
filter_list2 = np.ndarray.tolist(cell_filter2.reshape(1,-1)[0])
filter_list = [x - 1 for x in filter_list]
filter_list2 = [x - 1 for x in filter_list2]
filter_list = [elem for elem in filter_list if elem > 0]
filter_list2 = [elem for elem in filter_list2 if elem > 0]
mask_data2 = data2.get('mask_track')



#relevant cell shape and motility parameter values
v =np.vstack((mask_data[filter_list,:,5],mask_data2[filter_list2,:,5]))
#v = np.nan_to_num(v)
area = np.vstack((mask_data[filter_list,:,12],mask_data2[filter_list2,:,12]))
perim = np.vstack((mask_data[filter_list,:,17] , mask_data2[filter_list2,:,17]))
ecc =np.vstack((mask_data[filter_list,:,16] ,mask_data2[filter_list2,:,16]))
circ = np.vstack((mask_data[filter_list,:,15] , mask_data2[filter_list2,:,15]))

mean_v = np.nanmean(v,1)
top_v = np.percentile(mean_v, 50, axis=0)
fastCellIndex = np.where(mean_v>top_v)[0]
#v = v[fastCellIndex]
#area = area[fastCellIndex]
#perim = perim[fastCellIndex]
#ecc = ecc[fastCellIndex]
#circ = circ[fastCellIndex]

cell_shape = perim / np.sqrt(area)
mean_ecc = np.nanmean(ecc,1)
mean_circ = np.nanmean(circ,1)
coord_front = mask_data[:,:,6]
coord_lateral = mask_data[:,:,7]
coord_back = mask_data[:,:,8]
neighbours = mask_data[:,:,11]
coord_nb = mask_data[:,:,9]
adj_data = scp.io.loadmat('/mnt/d/Jessica/adjacency_flattened_3_8_1.mat')
adj_mat1 = adj_data.get('adj_mat')
adj_mat1 = adj_mat1[np.ix_(filter_list,filter_list)]
adj_data2 = scp.io.loadmat('/mnt/d/Jessica/adjacency_flattened_3_8_2.mat')
adj_mat2 = adj_data2.get('adj_mat')
adj_mat2 = adj_mat2[np.ix_(filter_list2,filter_list2)]

adj_mat = np.empty((adj_mat1.shape[0]+adj_mat2.shape[0],adj_mat1.shape[0]+adj_mat2.shape[0]))
adj_mat[:] = np.nan
adj_mat[:adj_mat1.shape[0],:adj_mat1.shape[0]] = adj_mat1
adj_mat[adj_mat1.shape[0]:adj_mat1.shape[0]+adj_mat2.shape[0],adj_mat1.shape[0]:adj_mat1.shape[0]+adj_mat2.shape[0]] =adj_mat2
random.seed(42)

nonan_arr = np.argwhere((~np.isnan(v[:,0:frame-2])).all(axis=1))
non_nan_shape = np.argwhere((~np.isnan(cell_shape[:,0:frame])).all(axis=1))

v_not_nan = ~np.isnan(v[:,0:178]).any(axis=1)
area_not_nan = ~np.isnan(area).any(axis=1)
#ecc_not_nan = ~np.isnan(ecc).any(axis=1)
circ_not_nan = ~np.isnan(circ).any(axis=1)
perim_not_nan = ~np.isnan(perim).any(axis=1)
intersect = np.logical_and(v_not_nan, area_not_nan)
#intersect = np.logical_and(intersect, ecc_not_nan)
intersect = np.logical_and(intersect, circ_not_nan)
intersect = np.logical_and(intersect, perim_not_nan)
inter_idx = np.where(intersect)[0]

print(inter_idx)
sample_cells = inter_idx.tolist()#sample(inter_idx.tolist(),280)
adj_mat = adj_mat[np.ix_(sample_cells,sample_cells)]
v_sampled = (v[sample_cells,0:frame-1]).T
shape_sampled = (cell_shape[sample_cells,0:frame-1]).T
area_sampled = (area[sample_cells,0:frame-1]).T
circ_sampled = (circ[sample_cells,0:frame-1]).T
samples = np.dstack((v_sampled,shape_sampled,area_sampled,circ_sampled))
#samples = np.dstack((v_sampled,shape_sampled,area_sampled))
#samples = np.asarray(v_sampled)
plt.figure(figsize=(18, 6))
plt.plot(v_sampled[:,[0, 1]])
plt.legend(["cell_0", "cell_1"])
plt.show()

plt.figure(figsize=(8, 8))
df_sample = pd.DataFrame(v_sampled.T)
sns.heatmap(df_sample.corr())
#plt.matshow(df_sample.corr(), 0)
plt.xlabel("cell_number")
plt.ylabel("cell_number")
plt.show()

df_shape = pd.DataFrame(shape_sampled.T)
sns.heatmap(df_shape.corr())
#plt.matshow(df_sample.corr(), 0)
plt.xlabel("cell_number")
plt.ylabel("cell_number")
plt.show()

train_size, val_size = 0.6, 0.2


def preprocess(data_array: np.ndarray, train_size: float, val_size: float):

    num_time_steps = data_array.shape[0]
    num_train, num_val = (
        int(num_time_steps * train_size),
        int(num_time_steps * val_size),
    )
    train_array = data_array[:num_train]
    #train_v = train_array
    train_v = train_array[:,:,0]
    train_shape = train_array[:,:,1]
    train_area = train_array[:,:,2]
    train_circ = train_array[:,:,3]

    mean, std = train_shape.mean(axis=0), train_shape.std(axis=0)
    mean_a, std_a = train_area.mean(axis=0), train_area.std(axis=0)
    mean_c, std_c = train_circ.mean(axis=0), train_circ.std(axis=0)

    train_shape = (train_shape - mean) / std
    val_shape= (data_array[num_train : (num_train + num_val),:,1] - mean) / std
    test_shape = (data_array[(num_train + num_val) :,:,1] - mean) / std

    train_area = (train_area - mean_a) / std_a
    val_area = (data_array[num_train: (num_train + num_val), :, 2] - mean_a) / std_a
    test_area = (data_array[(num_train + num_val):, :, 2] - mean_a) / std_a

    train_circ = (train_circ - mean_c) / std_c
    val_circ = (data_array[num_train: (num_train + num_val), :, 3] - mean_c) / std_c
    test_circ = (data_array[(num_train + num_val):, :, 3] - mean_c) / std_c

    val_v = data_array[num_train: (num_train + num_val), :,0]
    test_v = data_array[(num_train + num_val):, :,0]
    #val_v = (data_array[num_train: (num_train + num_val), :, 0] - mean_v) / std_v
    #test_v = (data_array[(num_train + num_val):, :, 0] - mean_v) / std_v
    train_array[:,:,0] = train_v
    train_array[:,:,1] = train_shape
    train_array[:,:,2] = train_area
    train_array[:, :, 3] = train_circ
    #val_array = np.dstack((val_v,val_shape))#,val_area,val_circ))
    #test_array = np.dstack((test_v,test_shape))#,test_area,test_circ))
    val_array = np.dstack((val_v,val_shape,val_area,val_circ))
    test_array = np.dstack((test_v,test_shape,test_area,test_circ))
    return train_array, val_array, test_array


train_array, val_array, test_array = preprocess(samples, train_size, val_size)
#trainX, valX,testX = preprocess(shape_sampled,train_size,val_size)
print(f"train set size: {train_array.shape}")
print(f"validation set size: {val_array.shape}")
print(f"test set size: {test_array.shape}")

batch_size = 4
input_sequence_length = 6
forecast_horizon = 3
multi_horizon = False


def create_tf_dataset(
    data_array: np.ndarray,
    input_sequence_length: int,
    forecast_horizon: int,
    batch_size: int = 4,
    shuffle=True,
    multi_horizon=False,
):


    inputs = keras.utils.timeseries_dataset_from_array(
        data_array[:-forecast_horizon],
        None,
        sequence_length=input_sequence_length,
        shuffle=False,
        batch_size=batch_size,
    )

    target_offset = (
        input_sequence_length
        if multi_horizon
        else input_sequence_length + forecast_horizon - 1
    )
    target_seq_length = forecast_horizon if multi_horizon else 1
    targets = keras.utils.timeseries_dataset_from_array(
        data_array[target_offset:,:,0],
        None,
        sequence_length=target_seq_length,
        shuffle=False,
        batch_size=batch_size,
    )

    dataset = tf.data.Dataset.zip((inputs, targets))
    if shuffle:
        dataset = dataset.shuffle(100)

    return dataset.prefetch(16).cache()


train_dataset, val_dataset = (
    create_tf_dataset(data_array, input_sequence_length, forecast_horizon, batch_size)
    for data_array in [train_array, val_array]
)

test_dataset = create_tf_dataset(
    test_array,
    input_sequence_length,
    forecast_horizon,
    batch_size=test_array.shape[0],
    shuffle=False,
    multi_horizon=multi_horizon,
)




def compute_adjacency_matrix(
    adj: np.ndarray
):

    mask = np.isfinite(adj)
    mask[mask] = adj[mask] > 0
    return mask


class GraphInfo:
    def __init__(self, edges: typing.Tuple[list, list], num_nodes: int):
        self.edges = edges
        self.num_nodes = num_nodes


adjacency_matrix = compute_adjacency_matrix(adj_mat)



class GraphConv(layers.Layer):
    def __init__(
        self,
        in_feat,
        out_feat,
        graph_info: GraphInfo,
        aggregation_type="mean",
        combination_type="concat",
        activation: typing.Optional[str] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.in_feat = in_feat
        self.out_feat = out_feat
        self.graph_info = graph_info
        self.aggregation_type = aggregation_type
        self.combination_type = combination_type
        self.weight = self.add_weight(
            initializer=keras.initializers.GlorotUniform(),
            shape=(in_feat, out_feat),
            dtype="float32",
            trainable=True,
        )
        self.activation = layers.Activation(activation)

    def aggregate(self, neighbour_representations):
        aggregation_func = {
            "sum": tf.math.unsorted_segment_sum,
            "mean": tf.math.unsorted_segment_mean,
            "max": tf.math.unsorted_segment_max,
        }.get(self.aggregation_type)

        if aggregation_func:
            return aggregation_func(
                neighbour_representations,
                self.graph_info.edges[0],
                num_segments=self.graph_info.num_nodes,
            )

        raise ValueError(f"Invalid aggregation type: {self.aggregation_type}")

    def compute_nodes_representation(self, features):
        return ops.matmul(features, self.weight)

    def compute_aggregated_messages(self, features):
        neighbour_representations = tf.gather(features, self.graph_info.edges[1])
        aggregated_messages = self.aggregate(neighbour_representations)
        return ops.matmul(aggregated_messages, self.weight)

    def update(self, nodes_representation, aggregated_messages):
        if self.combination_type == "concat":
            h = ops.concatenate([nodes_representation, aggregated_messages], axis=-1)
        elif self.combination_type == "add":
            h = nodes_representation + aggregated_messages
        else:
            raise ValueError(f"Invalid combination type: {self.combination_type}.")
        return self.activation(h)

    def call(self, features):
        nodes_representation = self.compute_nodes_representation(features)
        aggregated_messages = self.compute_aggregated_messages(features)
        return self.update(nodes_representation, aggregated_messages)

class LSTMGC(layers.Layer):
    """Layer comprising a convolution layer followed by LSTM and dense layers."""

    def __init__(
        self,
        in_feat,
        out_feat,
        lstm_units: int,
        input_seq_len: int,
        output_seq_len: int,
        adjacency_matrix,
        graph_conv_params: typing.Optional[dict] = None,
        **kwargs,
    ):
        super().__init__(**kwargs)

        # graph conv layer
        if graph_conv_params is None:
            graph_conv_params = {
                "aggregation_type": "mean",
                "combination_type": "concat",
                "activation": None,
            }

        node_indices, neighbor_indices = np.where(adjacency_matrix == 1)
        graph_info = GraphInfo(
            edges=(node_indices.tolist(), neighbor_indices.tolist()),
            num_nodes=adjacency_matrix.shape[0],
        )
        print(f"number of nodes: {graph_info.num_nodes}, number of edges: {len(graph_info.edges[0])}")
        self.graph_conv = GraphConv(in_feat, out_feat, graph_info, **graph_conv_params)

        self.lstm = layers.LSTM(lstm_units, activation="tanh")

        self.dense = layers.Dense(output_seq_len)

        self.input_seq_len, self.output_seq_len = input_seq_len, output_seq_len

    def call(self, inputs):
        # convert shape to  (num_nodes, batch_size, input_seq_len, in_feat)
        inputs = ops.transpose(inputs, [2, 0, 1, 3])

        gcn_out = self.graph_conv(
            inputs
        )  # gcn_out has shape: (num_nodes, batch_size, input_seq_len, out_feat)
        shape = ops.shape(gcn_out)
        num_nodes, batch_size, input_seq_len, out_feat = (
            shape[0],
            shape[1],
            shape[2],
            shape[3],
        )

        # LSTM takes only 3D tensors as input
        gcn_out = ops.reshape(
            gcn_out, (batch_size * num_nodes, input_seq_len, out_feat)
        )
        lstm_out = self.lstm(
            gcn_out
        )  # lstm_out has shape: (batch_size * num_nodes, lstm_units)

        dense_output = self.dense(
            lstm_out
        )  # dense_output has shape: (batch_size * num_nodes, output_seq_len)
        output = ops.reshape(dense_output, (num_nodes, batch_size, self.output_seq_len))
        return ops.transpose(
            output, [1, 2, 0]
        )  # returns Tensor of shape (batch_size, output_seq_len, num_nodes)


"""
## Model training
"""

#in_feat = 4
in_feat = 4
batch_size = 4
epochs = 100
input_sequence_length = 6
forecast_horizon = 3
multi_horizon = False
out_feat = 10
lstm_units = 128
graph_conv_params = {
    "aggregation_type": "mean",
    "combination_type": "concat",
    "activation": 'relu',
}

st_gcn = LSTMGC(
    in_feat,
    out_feat,
    lstm_units,
    input_sequence_length,
    forecast_horizon,
    adjacency_matrix,
    graph_conv_params,
)
inputs = layers.Input((input_sequence_length, adjacency_matrix.shape[0], in_feat))
outputs = st_gcn(inputs)

model = keras.models.Model(inputs, outputs)

model.compile(
    optimizer=keras.optimizers.RMSprop(learning_rate=0.0007),
    loss=keras.losses.MeanSquaredError()
)
history = model.fit(
    train_dataset,
    validation_data=val_dataset,
    epochs=epochs,
    callbacks=[keras.callbacks.EarlyStopping(patience=15)]
)
print(history.history.keys())

# summarize history for loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['Training Loss', 'Validation Loss'], loc='upper right')
plt.show()
#model.save('gcn_lstm_model.keras')
#

x_test, y = next(test_dataset.as_numpy_iterator())
y_pred = model.predict(x_test,verbose=1)
#y_pred_bool = np.argmax(y_pred, axis=1)
plt.figure(figsize=(18, 6))
plt.plot(y[:, 0, 0])
plt.plot(y_pred[:, 0, 0])
plt.legend(["actual", "predicted"])
plt.xlabel('Frame')
plt.ylabel('Velocity')


model_mse = np.square(y_pred[:, 0, :] - y[:, 0, :]).mean()
#print(f"naive MSE: {naive_mse},
print(f"model MSE: {model_mse}")


plt.figure(figsize=(18, 6))
plt.plot(y[:,0,62],label='Observed')
plt.plot(y_pred[:,0,62],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()

plt.figure(figsize=(18, 6))
plt.plot(y[:,0,6],label='Observed')
plt.plot(y_pred[:,0,6],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()


plt.figure(figsize=(18, 6))
plt.plot(y[:,0,8],label='Observed')
plt.plot(y_pred[:,0,8],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()

plt.figure(figsize=(18, 6))
plt.plot(y[:,0,18],label='Observed')
plt.plot(y_pred[:,0,18],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()

plt.figure(figsize=(18, 6))
plt.plot(y[:,0,19],label='Observed')
plt.plot(y_pred[:,0,19],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()

plt.figure(figsize=(18, 6))
plt.plot(y[:,0,17],label='Observed')
plt.plot(y_pred[:,0,17],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()

plt.figure(figsize=(18, 6))
plt.plot(y[:,0,12],label='Observed')
plt.plot(y_pred[:,0,12],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()

plt.figure(figsize=(18, 6))
plt.plot(y[:,0,15],label='Observed')
plt.plot(y_pred[:,0,15],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()

plt.figure(figsize=(18, 6))
plt.plot(y[:,0,32],label='Observed')
plt.plot(y_pred[:,0,32],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()


plt.figure(figsize=(18, 6))
plt.plot(y[:,0,45],label='Observed')
plt.plot(y_pred[:,0,45],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()


plt.figure(figsize=(18, 6))
plt.plot(y[:,0,27],label='Observed')
plt.plot(y_pred[:,0,27],label='Predicted')
plt.legend()
plt.xlabel('Frame')
plt.ylabel('Velocity')
plt.show()



