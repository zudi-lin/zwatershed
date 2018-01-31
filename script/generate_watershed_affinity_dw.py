# -*- coding: utf-8 -*-
#!/usr/bin/env python
import sys
import zwatershed
import h5py
import numpy as np

def readh5(filename, datasetname):
    fid=h5py.File(filename)
    # avoid wrap around: 1.1 -> 0 
    if datasetname not in fid.keys():
        print 'warning:',datasetname+' not found..'
        datasetname=fid.keys()[0]
    # v0. float
    data=np.array(fid[datasetname]).astype(np.float32)
    fid.close()
    return data
    
if __name__ == "__main__" :
  
    prediction_file=sys.argv[1]
    prediction_dataset= sys.argv[2]
    watershed_file= sys.argv[3]
    thd = [int(x) for x in sys.argv[4].split('_')]	
    aff_thd = [float(x) for x in sys.argv[5].split('_')]	
    dust_thd = int(sys.argv[6])
    
    #mask_file = prediction_file[:-7]+'mask.h5'
    p = readh5(prediction_file, prediction_dataset)
    print "Performing watershed"
    zwatershed.zwatershed_dw(p,thd, aff_thd, dust_thd, seg_save_path=watershed_file)
