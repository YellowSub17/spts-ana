import h5py 
import numpy as np 

from .filtering import ComboRun_Filters
from .plotting import ComboRun_Plots

import copy

class ComboRun(ComboRun_Filters, ComboRun_Plots):

    def __init__(self, fnames, hit_type='all'):

        
        self.fnames = fnames
        self.hit_type = hit_type #single, multi, all

        self.peak_xs = np.array([], dtype=int)
        self.peak_ys = np.array([], dtype=int)
        self.peak_is = np.array([], dtype=int)

        self.fname_inds = np.array([], dtype=int)
        self.img_inds = np.array([], dtype=int)
        
        
        for i, fname in enumerate(self.fnames):
            
        
            with h5py.File(fname, 'r')  as f:
                
                #number of peaks in each frame. almost all 0s, size about 2000
                run_npeaks = f['5_detect/n'][:]
                run_hits_loc =  np.where(run_npeaks>0)[0] #both single and multi
                #indices in the h5 that match the peak
                run_img_inds = np.array([file_i for file_i, npeaks_in_file_i in zip(run_hits_loc, run_npeaks[run_npeaks>0]) for _ in range(npeaks_in_file_i)])

                if len(run_hits_loc) == 0:
                    print(f'Warning: file {fname} has no hits')
                
                #x, y, and peak intensities
                # full array is size about 2000 (frames) x 2 (at most peaks detected per frame).
                #ravel to make a 1d array
                run_peak_xs = f['5_detect/x'][:].ravel()
                run_peak_ys = f['5_detect/y'][:].ravel()
                run_peak_is = f['6_analyse/peak_sum'][:].ravel()

            
            #in frames where there is no peaks, the xy positions are -1.
            xygt0_loc = np.logical_and(run_peak_xs>0, run_peak_ys>0)
        
            #remove the non-peaks
            run_peak_is = run_peak_is[xygt0_loc]
            run_peak_xs = run_peak_xs[xygt0_loc]
            run_peak_ys = run_peak_ys[xygt0_loc]

    
            # for some reason, this filter needs to be don seperate to the one above
            igt0 = run_peak_is>0
            
            #remove the non-peaks
            run_peak_is = run_peak_is[igt0]
            run_peak_xs = run_peak_xs[igt0]
            run_peak_ys = run_peak_ys[igt0]
            run_img_inds = run_img_inds[igt0]
        
            #indices that match to self.fname_inds
            run_fname_inds = np.zeros(run_peak_is.size, dtype=int)+i
          
            self.fname_inds = np.concatenate((run_fname_inds, self.fname_inds))
            self.img_inds = np.concatenate((run_img_inds, self.img_inds))
            self.peak_xs = np.concatenate((run_peak_xs, self.peak_xs))
            self.peak_ys = np.concatenate((run_peak_ys, self.peak_ys))
            self.peak_is = np.concatenate((run_peak_is, self.peak_is))


        #initialize running filter
        self.filter = np.ones(len(self.peak_is)).astype(bool)


    def copy(self):
        return copy.deepcopy(self)


def generate_filenames(ranges, a, let='a'):
    files = []
    for r in ranges:
        files +=[f'/home/tejvarmay/scattering_data/data/newdata/data00{i}_analysis_{let}{a}/spts.cxi' for i in range(r[0], r[1])]
    return files
