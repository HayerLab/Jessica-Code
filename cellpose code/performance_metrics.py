# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 11:59:06 2023

Code modified for use from:
https://github.com/HenriquesLab/ZeroCostDL4Mic/blob/master/Colab_notebooks/Cellpose_2D_ZeroCostDL4Mic.ipynb
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import cellpose
from cellpose import models
from numba import jit
from tqdm import tqdm
from scipy.optimize import linear_sum_assignment
from collections import namedtuple
from tifffile import imwrite, imread
import pandas as pd
import csv
from glob import glob
from astropy.visualization import simple_norm
import torch

# model name and path
#@markdown ###Do you want to assess the model you just trained ?
Use_the_current_trained_model = True #@param {type:"boolean"}

#@markdown ###If not, indicate which model you want to assess:

QC_model = "Own_model" #@param ["Cytoplasm","Cytoplasm2","LiveCell", "TissueNet", "Cytoplasm2_Omnipose", "Bacteria_Omnipose", "Nuclei", "Own_model"]

#@markdown ###If using your own model, please provide the path to the model (not the folder):

QC_model_path = "E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig tiff\\training\\models\\Improved_HUVEC_monolayer" #@param {type:"string"}

#@markdown ###If using the Cellpose or Omnipose models, please indicate where you want to save the results:
Saving_path = "E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig tiff\\test" #@param {type:"string"}
QC_model_folder = "E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig tiff\\training\\models"



# Here we load the def that perform the QC, code taken from StarDist  https://github.com/mpicbg-csbd/stardist/blob/master/stardist/matching.py



matching_criteria = dict()

def label_are_sequential(y):
    """ returns true if y has only sequential labels from 1... """
    labels = np.unique(y)
    return (set(labels)-{0}) == set(range(1,1+labels.max()))


def is_array_of_integers(y):
    return isinstance(y,np.ndarray) and np.issubdtype(y.dtype, np.integer)


def _check_label_array(y, name=None, check_sequential=False):
    err = ValueError("{label} must be an array of {integers}.".format(
        label = 'labels' if name is None else name,
        integers = ('sequential ' if check_sequential else '') + 'non-negative integers',
    ))
    is_array_of_integers(y) or print("An error occured")
    if check_sequential:
        label_are_sequential(y) or print("An error occured")
    else:
        y.min() >= 0 or print("An error occured")
    return True


def label_overlap(x, y, check=True):
    if check:
        _check_label_array(x,'x',True)
        _check_label_array(y,'y',True)
        x.shape == y.shape #or (raise ValueError("x and y must have the same shape"))
    return _label_overlap(x, y)

@jit(nopython=True)
def _label_overlap(x, y):
    x = x.ravel()
    y = y.ravel()
    overlap = np.zeros((1+x.max(),1+y.max()), dtype=np.uint)
    for i in range(len(x)):
        overlap[x[i],y[i]] += 1
    return overlap


def intersection_over_union(overlap):
    _check_label_array(overlap,'overlap')
    if np.sum(overlap) == 0:
        return overlap
    n_pixels_pred = np.sum(overlap, axis=0, keepdims=True)
    n_pixels_true = np.sum(overlap, axis=1, keepdims=True)
    return overlap / (n_pixels_pred + n_pixels_true - overlap)

matching_criteria['iou'] = intersection_over_union


def intersection_over_true(overlap):
    _check_label_array(overlap,'overlap')
    if np.sum(overlap) == 0:
        return overlap
    n_pixels_true = np.sum(overlap, axis=1, keepdims=True)
    return overlap / n_pixels_true

matching_criteria['iot'] = intersection_over_true


def intersection_over_pred(overlap):
    _check_label_array(overlap,'overlap')
    if np.sum(overlap) == 0:
        return overlap
    n_pixels_pred = np.sum(overlap, axis=0, keepdims=True)
    return overlap / n_pixels_pred

matching_criteria['iop'] = intersection_over_pred


def precision(tp,fp,fn):
    return tp/(tp+fp) if tp > 0 else 0
def recall(tp,fp,fn):
    return tp/(tp+fn) if tp > 0 else 0
def accuracy(tp,fp,fn):
    # also known as "average precision" (?)
    # -> https://www.kaggle.com/c/data-science-bowl-2018#evaluation
    return tp/(tp+fp+fn) if tp > 0 else 0
def f1(tp,fp,fn):
    # also known as "dice coefficient"
    return (2*tp)/(2*tp+fp+fn) if tp > 0 else 0


def _safe_divide(x,y):
    return x/y if y>0 else 0.0

def matching(y_true, y_pred, thresh=0.5, criterion='iou', report_matches=False):
    """Calculate detection/instance segmentation metrics between ground truth and predicted label images.
    Currently, the following metrics are implemented:
    'fp', 'tp', 'fn', 'precision', 'recall', 'accuracy', 'f1', 'criterion', 'thresh', 'n_true', 'n_pred', 'mean_true_score', 'mean_matched_score', 'panoptic_quality'
    Corresponding objects of y_true and y_pred are counted as true positives (tp), false positives (fp), and false negatives (fn)
    whether their intersection over union (IoU) >= thresh (for criterion='iou', which can be changed)
    * mean_matched_score is the mean IoUs of matched true positives
    * mean_true_score is the mean IoUs of matched true positives but normalized by the total number of GT objects
    * panoptic_quality defined as in Eq. 1 of Kirillov et al. "Panoptic Segmentation", CVPR 2019
    Parameters
    ----------
    y_true: ndarray
        ground truth label image (integer valued)
        predicted label image (integer valued)
    thresh: float
        threshold for matching criterion (default 0.5)
    criterion: string
        matching criterion (default IoU)
    report_matches: bool
        if True, additionally calculate matched_pairs and matched_scores (note, that this returns even gt-pred pairs whose scores are below  'thresh')
    Returns
    """
    _check_label_array(y_true,'y_true')
    _check_label_array(y_pred,'y_pred')
    y_true.shape == y_pred.shape #or _raise(ValueError("y_true ({y_true.shape}) and y_pred ({y_pred.shape}) have different shapes".format(y_true=y_true, y_pred=y_pred)))
    criterion in matching_criteria #or _raise(ValueError("Matching criterion '%s' not supported." % criterion))
    if thresh is None: thresh = 0
    thresh = float(thresh) if np.isscalar(thresh) else map(float,thresh)

    y_true, _, map_rev_true = relabel_sequential(y_true)
    y_pred, _, map_rev_pred = relabel_sequential(y_pred)

    overlap = label_overlap(y_true, y_pred, check=False)
    scores = matching_criteria[criterion](overlap)
    assert 0 <= np.min(scores) <= np.max(scores) <= 1

    # ignoring background
    scores = scores[1:,1:]
    n_true, n_pred = scores.shape
    n_matched = min(n_true, n_pred)

    def _single(thr):
        not_trivial = n_matched > 0 and np.any(scores >= thr)
        if not_trivial:
            # compute optimal matching with scores as tie-breaker
            costs = -(scores >= thr).astype(float) - scores / (2*n_matched)
            true_ind, pred_ind = linear_sum_assignment(costs)
            assert n_matched == len(true_ind) == len(pred_ind)
            match_ok = scores[true_ind,pred_ind] >= thr
            tp = np.count_nonzero(match_ok)
        else:
            tp = 0
        fp = n_pred - tp
        fn = n_true - tp

        # the score sum over all matched objects (tp)
        sum_matched_score = np.sum(scores[true_ind,pred_ind][match_ok]) if not_trivial else 0.0

        # the score average over all matched objects (tp)
        mean_matched_score = _safe_divide(sum_matched_score, tp)
        # the score average over all gt/true objects
        mean_true_score    = _safe_divide(sum_matched_score, n_true)
        panoptic_quality   = _safe_divide(sum_matched_score, tp+fp/2+fn/2)

        stats_dict = dict (
            criterion          = criterion,
            thresh             = thr,
            fp                 = fp,
            tp                 = tp,
            fn                 = fn,
            precision          = precision(tp,fp,fn),
            recall             = recall(tp,fp,fn),
            accuracy           = accuracy(tp,fp,fn),
            f1                 = f1(tp,fp,fn),
            n_true             = n_true,
            n_pred             = n_pred,
            mean_true_score    = mean_true_score,
            mean_matched_score = mean_matched_score,
            panoptic_quality   = panoptic_quality,
        )
        if bool(report_matches):
            if not_trivial:
                stats_dict.update (
                    # int() to be json serializable
                    matched_pairs  = tuple((int(map_rev_true[i]),int(map_rev_pred[j])) for i,j in zip(1+true_ind,1+pred_ind)),
                    matched_scores = tuple(scores[true_ind,pred_ind]),
                    matched_tps    = tuple(map(int,np.flatnonzero(match_ok))),
                )
            else:
                stats_dict.update (
                    matched_pairs  = (),
                    matched_scores = (),
                    matched_tps    = (),
                )
        return namedtuple('Matching',stats_dict.keys())(*stats_dict.values())

    return _single(thresh) if np.isscalar(thresh) else tuple(map(_single,thresh))



def matching_dataset(y_true, y_pred, thresh=0.5, criterion='iou', by_image=False, show_progress=True, parallel=False):
    """matching metrics for list of images, see `stardist.matching.matching`
    """
    len(y_true) == len(y_pred)
    return matching_dataset_lazy (
        tuple(zip(y_true,y_pred)), thresh=thresh, criterion=criterion, by_image=by_image, show_progress=show_progress, parallel=parallel,
    )



def matching_dataset_lazy(y_gen, thresh=0.5, criterion='iou', by_image=False, show_progress=True, parallel=False):

    expected_keys = set(('fp', 'tp', 'fn', 'precision', 'recall', 'accuracy', 'f1', 'criterion', 'thresh', 'n_true', 'n_pred', 'mean_true_score', 'mean_matched_score', 'panoptic_quality'))

    single_thresh = False
    if np.isscalar(thresh):
        single_thresh = True
        thresh = (thresh,)

    tqdm_kwargs = {}
    tqdm_kwargs['disable'] = not bool(show_progress)
    if int(show_progress) > 1:
        tqdm_kwargs['total'] = int(show_progress)

    # compute matching stats for every pair of label images
    if parallel:
        from concurrent.futures import ThreadPoolExecutor
        fn = lambda pair: matching(*pair, thresh=thresh, criterion=criterion, report_matches=False)
        with ThreadPoolExecutor() as pool:
            stats_all = tuple(pool.map(fn, tqdm(y_gen,**tqdm_kwargs)))
    else:
        stats_all = tuple (
            matching(y_t, y_p, thresh=thresh, criterion=criterion, report_matches=False)
            for y_t,y_p in tqdm(y_gen,**tqdm_kwargs)
        )

    # accumulate results over all images for each threshold separately
    n_images, n_threshs = len(stats_all), len(thresh)
    accumulate = [{} for _ in range(n_threshs)]
    for stats in stats_all:
        for i,s in enumerate(stats):
            acc = accumulate[i]
            for k,v in s._asdict().items():
                if k == 'mean_true_score' and not bool(by_image):
                    # convert mean_true_score to "sum_matched_score"
                    acc[k] = acc.setdefault(k,0) + v * s.n_true
                else:
                    try:
                        acc[k] = acc.setdefault(k,0) + v
                    except TypeError:
                        pass

    # normalize/compute 'precision', 'recall', 'accuracy', 'f1'
    for thr,acc in zip(thresh,accumulate):
        set(acc.keys()) == expected_keys #or _raise(ValueError("unexpected keys"))
        acc['criterion'] = criterion
        acc['thresh'] = thr
        acc['by_image'] = bool(by_image)
        if bool(by_image):
            for k in ('precision', 'recall', 'accuracy', 'f1', 'mean_true_score', 'mean_matched_score', 'panoptic_quality'):
                acc[k] /= n_images
        else:
            tp, fp, fn, n_true = acc['tp'], acc['fp'], acc['fn'], acc['n_true']
            sum_matched_score = acc['mean_true_score']

            mean_matched_score = _safe_divide(sum_matched_score, tp)
            mean_true_score    = _safe_divide(sum_matched_score, n_true)
            panoptic_quality   = _safe_divide(sum_matched_score, tp+fp/2+fn/2)

            acc.update(
                precision          = precision(tp,fp,fn),
                recall             = recall(tp,fp,fn),
                accuracy           = accuracy(tp,fp,fn),
                f1                 = f1(tp,fp,fn),
                mean_true_score    = mean_true_score,
                mean_matched_score = mean_matched_score,
                panoptic_quality   = panoptic_quality,
            )

    accumulate = tuple(namedtuple('DatasetMatching',acc.keys())(*acc.values()) for acc in accumulate)
    return accumulate[0] if single_thresh else accumulate



# copied from scikit-image master for now (remove when part of a release)
def relabel_sequential(label_field, offset=1):

    offset = int(offset)
    if offset <= 0:
        raise ValueError("Offset must be strictly positive.")
    if np.min(label_field) < 0:
        raise ValueError("Cannot relabel array that contains negative values.")
    max_label = int(label_field.max()) # Ensure max_label is an integer
    if not np.issubdtype(label_field.dtype, np.integer):
        new_type = np.min_scalar_type(max_label)
        label_field = label_field.astype(new_type)
    labels = np.unique(label_field)
    labels0 = labels[labels != 0]
    new_max_label = offset - 1 + len(labels0)
    new_labels0 = np.arange(offset, new_max_label + 1)
    output_type = label_field.dtype
    required_type = np.min_scalar_type(new_max_label)
    if np.dtype(required_type).itemsize > np.dtype(label_field.dtype).itemsize:
        output_type = required_type
    forward_map = np.zeros(max_label + 1, dtype=output_type)
    forward_map[labels0] = new_labels0
    inverse_map = np.zeros(new_max_label + 1, dtype=output_type)
    inverse_map[offset:] = labels0
    relabeled = forward_map[label_field]
    return relabeled, forward_map, inverse_map

Target_QC_folder = "E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig tiff\\test" #@param{type:"string"}
Prediction_folder = "E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig tiff\\results" #@param{type:"string"}

# Segmentation parameters:
Object_diameter =  105

Flow_threshold = 0.4 

mask_threshold=0 

if QC_model == "Own_model":
  channels=[2,1] #green channel for VE-cadherin stain and red channel for nuclear stain

model = models.CellposeModel(gpu=True, device=torch.device('cuda'),
                             pretrained_model=QC_model_path)

# Here we need to make predictions

#for name in glob(os.path.join(Prediction_folder,'*.tiff')):
  
print("Performing prediction on: "+ Prediction_folder)
#image = io.imread(name)
image = glob(os.path.join(Prediction_folder, '*.tiff'))

images = [cellpose.io.imread(f) for f in image]
if QC_model == "Own_model":
    masks, flows, styles = model.eval(images, diameter=Object_diameter, flow_threshold=Flow_threshold, cellprob_threshold=mask_threshold, channels=[2,1])
else:
    masks, flows, styles, diams = models.eval(images, diameter=Object_diameter, flow_threshold=Flow_threshold, cellprob_threshold=mask_threshold, channels=channels)

# os.chdir(Saving_path)
for (i, mask) in enumerate(masks): 
    imwrite(Saving_path+"\\predictions_"+str(i+1)+".tif", mask) 
 # imsave(str(short_name[0])+".tif", masks)

# Here we start testing the differences between GT and predicted masks
true_segmented = glob('E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig tiff\\test\\*.npy')
images_segment = [cellpose.io.imread(f) for f in true_segmented]
for i,img in enumerate(images_segment):
    np.savetxt(Saving_path+"\\true_"+str(i)+".csv", img, delimiter=",")
    imwrite(Saving_path+"\\images_"+str(i+1)+".tif", img) 
with open(Prediction_folder+"\\checks\\"+QC_model+".csv", "w", newline='') as file:
  writer = csv.writer(file, delimiter=",")
  writer.writerow(["image","Prediction v. GT Intersection over Union", "false positive", "true positive", "false negative", "precision", "recall", "accuracy", "f1 score", "n_true", "n_pred", "mean_true_score", "mean_matched_score", "panoptic_quality"])  

# define the images
  i = 0
  files = cellpose.io.get_image_files(Saving_path, ('_seg'))
  print(files)
  images = [cellpose.io.imread(f) for f in files]

  for n in range(13):
    if not os.path.isdir(images[i]):
      print('Running QC on: '+str(n))
      
      test_input = images[i]
      test_prediction = cellpose.io.imread(Saving_path+"\\predictions_"+str(n+1)+".tif")
     # test_ground_truth_image = io.imread(os.path.join(Target_QC_folder, n))
      test_ground_truth_image = cellpose.io.imread(Saving_path+"\\images_"+str(n+1)+".tif")
      # Calculate the matching (with IoU threshold `thresh`) and all metrics

      stats = matching(test_ground_truth_image, test_prediction, thresh=0.5)
      

      #Convert pixel values to 0 or 255
      test_prediction_0_to_255 = test_prediction
      test_prediction_0_to_255[test_prediction_0_to_255>0] = 255

      #Convert pixel values to 0 or 255
      test_ground_truth_0_to_255 = test_ground_truth_image
      test_ground_truth_0_to_255[test_ground_truth_0_to_255>0] = 255


      # Intersection over Union metric

      intersection = np.logical_and(test_ground_truth_0_to_255, test_prediction_0_to_255)
      union = np.logical_or(test_ground_truth_0_to_255, test_prediction_0_to_255)
      iou_score =  np.sum(intersection) / np.sum(union)
      writer.writerow([n, str(iou_score), str(stats.fp), str(stats.tp), str(stats.fn), str(stats.precision), str(stats.recall), str(stats.accuracy), str(stats.f1), str(stats.n_true), str(stats.n_pred), str(stats.mean_true_score), str(stats.mean_matched_score), str(stats.panoptic_quality)])
      i+=1
from tabulate import tabulate

df = pd.read_csv (Prediction_folder+"\\checks\\"+QC_model+".csv")
print(tabulate(df, headers='keys', tablefmt='psql'))

df.to_csv('E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig tiff\\test\\performance_updated.csv')
# ------------- For display ------------
print('--------------------------------------------------------------')

def show_QC_results(n_list):
    index = 0
    for i in n_list: 
      plt.figure(figsize=(25,5))
      source_image = imread("E:\\Jessica\\segmentation\\3_7_1\\Composites\\aligned\\export_fig tiff\\test\\composite_aligned_"+str(i)+".tiff")
    
      target_image = cellpose.io.imread(Saving_path+"\\images_"+str(index+1)+".tif")
      prediction = cellpose.io.imread(Saving_path+"\\predictions_"+str(index+1)+".tif")
    
      #stats = matching(prediction, target_image, thresh=0.5)
    
      target_image_mask = np.empty_like(target_image)
      target_image_mask[target_image > 0] = 255
      target_image_mask[target_image == 0] = 0
      
      prediction_mask = np.empty_like(prediction)
      prediction_mask[prediction > 0] = 255
      prediction_mask[prediction == 0] = 0
    
      intersection = np.logical_and(target_image_mask, prediction_mask)
      union = np.logical_or(target_image_mask, prediction_mask)
      iou_score =  np.sum(intersection) / np.sum(union)
    
      #norm = simple_norm(source_image, percent = 99)
    
      #Input
      plt.subplot(1,4,1)
      plt.axis('off')
      #if n_channel > 1:
      plt.imshow(source_image)
      plt.title('Input')
    
      #Ground-truth
      plt.subplot(1,4,2)
      plt.axis('off')
      plt.imshow(target_image_mask, aspect='equal', cmap='Greens')
      plt.title('Ground Truth')
    
      #Prediction
      plt.subplot(1,4,3)
      plt.axis('off')
      plt.imshow(prediction_mask, aspect='equal', cmap='Purples')
      plt.title('Prediction')
    
      #Overlay
      plt.subplot(1,4,4)
      plt.axis('off')
      plt.imshow(target_image_mask, cmap='Greens')
      plt.imshow(prediction_mask, alpha=0.5, cmap='Purples')
      plt.title('Ground Truth and Prediction, Intersection over Union:'+str(round(iou_score,3 )));
      plt.savefig(Saving_path+'\\performance_'+str(i)+'.png',bbox_inches='tight',pad_inches=0)
      index += 1


show_QC_results([46,47,48,49,50,51,52,53,54,55,56,57,59])