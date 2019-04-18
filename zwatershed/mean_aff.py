import numpy as np
from scipy.sparse import coo_matrix
from scipy.ndimage import grey_erosion, grey_dilation

# mean agglomeration
SIX_CONNECTED = np.array([[[False, False, False],
                           [False, True, False],
                           [False, False, False]],
                          [[False, True, False],
                           [True, True, True],
                           [False, True, False]],
                          [[False, False, False],
                           [False, True, False],
                           [False, False, False]]])

def mean_agglomeration(seg, aff, thd):
    rg_ma = ma_get_region_graph(seg, aff)
    new_seg = ma_merge_region(seg, rg_ma, thd)
    return new_seg

def ma_get_region_graph(seg, aff):
    # input:
    # output: N*3, (seg_id_1, seg_id_2, mean_aff)
    mpred = aff.mean(axis=0)
    # find boundary 
    ws_eroded = grey_erosion(seg, footprint=SIX_CONNECTED)
    ws_dilated = grey_dilation(seg, footprint=SIX_CONNECTED)
    different = (ws_eroded != 0) & (ws_eroded != ws_dilated)
    id1 = ws_eroded[different]
    id2 = ws_dilated[different]
    id1id2pred = mpred[different]
    m = coo_matrix((id1id2pred, (id1, id2)))
    m.sum_duplicates()
    ma = coo_matrix((np.ones(len(id1)), (id1, id2)))
    ma.sum_duplicates()
    id1a, id2a = m.nonzero()
    mm = m.tocsr()
    mma = ma.tocsr()
    scores = mm[id1a, id2a] / mma[id1a, id2a]
    scores = scores.A1
    order = np.argsort(np.max(scores) - scores)
    scores = scores[order]
    # already sorted: id1a < id2a
    id1a = id1a[order]
    id2a = id2a[order]
    rg_ma = np.vstack((id1a,id2a,scores)).transpose(1,0)
    return rg_ma

def ma_merge_region(seg, rg_ma, thd):
    id1 = rg_ma[:,0].astype(seg.dtype)
    id2 = rg_ma[:,1].astype(seg.dtype)

    tomerge = np.where(rg_ma[:,2]>=float(thd))[0]
    print("#to merge ", len(tomerge))

    seg_new = merge_label(seg, id1[tomerge], id2[tomerge])
    return seg_new

def merge_label(seg, id1, id2):
    # relabel id2 into id1 
    labels = np.unique(seg)
    m = int(labels.max())                                                                
    label_map = np.zeros(m+1, int)                                                     
    label_map[labels] = labels
    label_map[id2] = id1
    """
    for row in range(len(id1)):
        seg[seg==id2[row]] = id1[row]
    """
    return label_map[seg]
 
