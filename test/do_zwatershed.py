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
    aff_thd = [x for x in sys.argv[5].split('_')]
    T_rel = int(aff_thd[-1]) 
    dust_thd = int(sys.argv[6])
    T_merge = sys.argv[7]
    
    print "Initial watershed"
    save_init = save_prefix+'_%s_%s_%s'%(aff_thd[0], aff_thd[1], aff_thd[3])+'.h5'
    save_out_pref = save_init[:-3]+'_%s_%d_%s'%(aff_thd[2],dust_thd,T_merge)
    save_out = [save_out_pref+'_'+str(x)+'.h5' for x in thd]
    if any([not os.path.exists(x) for x in save_out]):
        if not os.path.exists(save_init):
            print "compute zwater-init"
            p = readh5(prediction_file, prediction_dataset)
            aff_thres = map(float, aff_thd[:2])
            aff_thres = np.percentile(p, aff_thres) if T_rel else aff_thres
            out = zwatershed.zw_initial(p, aff_thres[0], aff_thres[1])
            p_prctile = np.percentile(p, np.arange(1,20)*0.05)
            writeh5(save_init,['seg','rg','counts','prctile'],
                    [out['seg'], out['rg'], out['counts'], p_prctile])
        else:
            print "load zwater-init"
            out=h5py.File(save_init)
            p_prctile = np.array(out['prctile'])
        seg = np.array(out['seg'])
        rg = np.array(out['rg']).reshape(-1,3)
        counts = np.array(out['counts'])
        T_aff_merge = p_prctile[int(np.floor(aff_thd[2]/0.05))] if T_rel else aff_thd[2]
        print "Performing watershed"
        zwatershed.zw_merge_region(seg, counts, rg,
                                      thd, T_aff_merge=0.2,T_dust=dust_thd,
                                      T_merge=float(T_merge),seg_save_path=save_out_pref)
