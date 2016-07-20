import numpy as np
import sys
import time
import os
import h5py
import os.path as op
import matplotlib.cm as cm
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
import array

sys.path.append('..')
from zwatershed import *
from edgelist_methods import *


def init():
    # -------------------------------- parameters ---------------------------------------
    path_to_folder = '/Users/chandansingh/drive/janelia/conv_net_scripts/'
    path_to_data = path_to_folder + 'data/'
    global segs_old, segs_new, rand_old, rand_new
    segs_old, segs_new, rand_old, rand_new = [[]], [[]], -1, -1
    threshes = [10, 2000]
    hdf5_gt_file = path_to_data + 'groundtruth_seg_thick.h5'  # /groups/turaga/home/turagas/data/FlyEM/fibsem_medulla_7col/tstvol-520-1-h5/groundtruth_seg_thick.h5'
    hdf5_pred_file = path_to_data + 'tstvol-1_2.h5'  # /tier2/turaga/singhc/train/output_200000/tstvol-1_2.h5'
    seg_save_path = path_to_data + 'out/'  # '/groups/turaga/home/singhc/evaluation/out/'
    seg_save_path_arb = path_to_data + 'out_arb/'  # '/groups/turaga/home/singhc/evaluation/out/'
    save_threshes = threshes
    p1, p2, p3 = 160, 170, 180  # 215, 214, 214 # 200, 200, 200

    # ----------------------------- load/shape data ------------------------------------
    # hdf5_gt = h5py.File(hdf5_gt_file, 'r')
    # hdf5_aff = h5py.File(hdf5_pred_file, 'r')
    # gt = np.asarray(hdf5_gt[hdf5_gt.keys()[0]], dtype='uint64')
    # aff = np.asarray(hdf5_aff[hdf5_aff.keys()[0]], dtype='float32')
    # aff = aff[:, p1:(-1 * p1), p2:(-1 * p2), p3:(-1 * p3)]
    # gt = trim_arbitrary_aff(gt, aff)
    gt = np.random.randint(0,high=30,size=(30,30,30)).astype('uint64')
    aff = np.random.rand(3,30,30,30).astype('float32')
    nhood = mknhood3d(1)
    node1, node2, edge_affs = affgraph_to_edgelist(aff, nhood)
    return (gt, aff, threshes, save_threshes, node1, node2, edge_affs, seg_save_path)


# ------------------------------ run tests -------------------------------------
def main():
    args = init()
    test_eval(args)
    # test_no_eval(args)
    # test_h5_eval(args)
    # test_h5_no_eval(args)


def print_final(segs_old, rand_old, segs_new, rand_new):
    print "--------Final--------"
    print rand_old
    print rand_new
    print "nsegs", len(np.unique(segs_old[0])), len(np.unique(segs_old[-1]))
    print "nsegs", len(np.unique(segs_old[0])), len(np.unique(segs_new[-1]))


# ------------------------------ test definitions -------------------------------------
def test_eval(args):
    (gt, aff, threshes, save_threshes, node1, node2, edge_affs, seg_save_path) = args
    print "\noriginal watershed..."
    start = time.clock()
    segs_old, rand_old = zwatershed_and_metrics(gt, aff, threshes, save_threshes)
    print "time: ", time.clock() - start

    print "\nnew watershed..."
    start = time.clock()
    segs_new, rand_new = zwatershed_and_metrics_arb(gt, node1, node2, edge_affs, threshes, save_threshes)
    print "time: ", time.clock() - start, "\n"
    return segs_old, rand_old, segs_new, rand_new


def test_no_eval(args):
    (gt, aff, threshes, save_threshes, node1, node2, edge_affs, seg_save_path) = args
    print "\noriginal watershed..."
    start = time.clock()
    segs_old = zwatershed(aff, threshes)
    print "time: ", time.clock() - start

    print "\nnew watershed..."
    start = time.clock()
    segs_new = zwatershed_arb(gt.shape, node1, node2, edge_affs, save_threshes)
    print "time: ", time.clock() - start, "\n"
    return segs_old, segs_new


def test_h5_eval(args):
    (gt, aff, threshes, save_threshes, node1, node2, edge_affs, seg_save_path) = args
    print "\noriginal watershed..."
    start = time.clock()
    rand_old = zwatershed_and_metrics_h5(gt, aff, threshes, save_threshes, seg_save_path)
    print "time: ", time.clock() - start

    print "\nnew watershed..."
    start = time.clock()
    rand_new = zwatershed_and_metrics_h5_arb(gt, node1, node2, edge_affs, threshes, save_threshes, seg_save_path_arb)
    print "time: ", time.clock() - start, "\n"
    return rand_old, rand_new


def test_h5_no_eval(args):
    (gt, aff, threshes, save_threshes, node1, node2, edge_affs, seg_save_path) = args
    print "\noriginal watershed..."
    start = time.clock()
    zwatershed_h5(aff, save_threshes, seg_save_path)
    print "time: ", time.clock() - start

    print "\nnew watershed..."
    start = time.clock()
    zwatershed_h5_arb(gt.shape, node1, node2, edge_affs, save_threshes, seg_save_path_arb)
    print "time: ", time.clock() - start, "\n"


if __name__ == '__main__':
    main()
