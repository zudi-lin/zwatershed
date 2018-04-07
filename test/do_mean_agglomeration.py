import numpy as np
from scipy.sparse import coo_matrix
from scipy.ndimage import grey_erosion, grey_dilation
import os,sys
import h5py
import zwatershed

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
            ma_rg = zwatershed.ma_get_region_graph(watershed, affinity)
            print "saving ", output_txt
            np.savetxt(output_txt, ma_rg, fmt='%.2f', delimiter=",")
        else:
            print "loading ", output_txt
            ma_rg = np.loadtxt(output_txt, delimiter=",")
        if watershed is None:
            watershed = readh5(sys.argv[3], sys.argv[4])
        new_seg = zwatershed.ma_merge_region(watershed, ma_rg, thd)
        writeh5(output_h5,'stack', new_seg)
