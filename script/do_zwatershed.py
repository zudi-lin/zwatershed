# -*- coding: utf-8 -*-
#!/usr/bin/env python
import os,sys
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

def writeh5(filename, datasetnames, dtarrays):
    fid=h5py.File(filename,'w')
    for k in range(len(datasetnames)):
        ds = fid.create_dataset(datasetnames[k], dtarrays[k].shape, compression="gzip", dtype=dtarrays[k].dtype)  
        ds[:] = dtarrays[k]
    fid.close()
    
if __name__ == "__main__" :
    prediction_file=sys.argv[1]
    prediction_dataset= sys.argv[2]
    save_prefix = sys.argv[3]
    thd = [int(x) for x in sys.argv[4].split('_')]	
    aff_thd = [float(x) for x in sys.argv[5].split('_')]	
    T_rel = int(aff_thd[-1]) 
    dust_thd = int(sys.argv[6])
    T_mst = int(sys.argv[7])
    T_merge = float(sys.argv[8])
    
    print "Initial watershed"
    save_init = save_prefix+str(aff_thd[0])+str(aff_thd[1])+str(aff_thd[3])+'.h5'
    p = None
    if not os.path.exists(save_init):
        p = readh5(prediction_file, prediction_dataset)
        aff_thres = np.percentile(p, aff_thd[:-1]) if T_rel else aff_thd
        out = zwatershed.zw_initial(p, aff_thres[0], aff_thres[1])
        p_prctile = np.percentile(p, np.arange(1,20)*0.05)
        writeh5(save_init,[],[out[],out[],out['counts'],p_prctile])
    print "Performing watershed"
    zwatershed.zwatershed(p,thd, aff_thd[:-1], T_aff_relative=T_rel,
                             T_dust=dust_thd, T_mst=T_mst,, T_merge=T_merge,
                             seg_save_path=watershed_file, 
                             )
