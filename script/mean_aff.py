import numpy as np
from scipy.sparse import coo_matrix
from scipy.ndimage import grey_erosion, grey_dilation
import os,sys
import h5py

def readh5(filename, datasetname='main'):
    fid=h5py.File(filename)
    # avoid wrap around: 1.1 -> 0 
    if datasetname not in fid.keys():
        print 'warning:',datasetname+' not found..'
        datasetname=fid.keys()[0]
    # v0. float
    data=np.array(fid[datasetname])
    fid.close()
    return data
    
def writeh5(filename, datasetname, dtarray):
    fid=h5py.File(filename,'w')
    ds = fid.create_dataset(datasetname, dtarray.shape, compression="gzip", dtype=dtarray.dtype)  
    ds[:] = dtarray 
    fid.close()
 

SIX_CONNECTED = np.array([[[False, False, False],
                           [False, True, False],
                           [False, False, False]],
                          [[False, True, False],
                           [True, True, True],
                           [False, True, False]],
                          [[False, False, False],
                           [False, True, False],
                           [False, False, False]]])
if __name__ == "__main__":
    # originial order: czyx
    output_txt = sys.argv[1]
    thd = sys.argv[2]
    output_h5 = output_txt[:-4]+'_'+thd+'.h5'

    if os.path.exists(output_h5):
        print "done already: ",output_h5
    else:
        watershed = None
        if not os.path.exists(output_txt):
            print "compute mean affinity"
            watershed = readh5(sys.argv[3], sys.argv[4])
            affinity = readh5(sys.argv[5], sys.argv[6]).astype(np.float32)
            mpred = affinity.mean(axis=0)
            # find boundary # find boundary # find boundary # find boundary # find boundary # find boundary # find boundary # find boundary # find boundary # find boundary # find boundary
            ws_eroded = grey_erosion(watershed, footprint=SIX_CONNECTED)
            ws_dilated = grey_dilation(watershed, footprint=SIX_CONNECTED)
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
            out = np.vstack((id1a,id2a,scores)).transpose(1,0)
            print "saving ", output_txt
            np.savetxt(output_txt, out, fmt='%.2f', delimiter=",")
        else:
            print "loading ", output_txt
            out = np.loadtxt(output_txt, delimiter=",")
        if watershed is None:
            watershed = readh5(sys.argv[3], sys.argv[4])
        outT = watershed.dtype
        tomerge = np.where(out[:,2]>=float(thd))[0]
        print "#to merge ", len(tomerge)
        for row in tomerge:
            watershed[watershed==out[row,1].astype(outT)] = out[row,0].astype(outT)
        print "saving ", output_h5
        writeh5(output_h5,'stack',watershed)
