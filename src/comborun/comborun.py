
from .filtering import ComboRun_Filters
from .plotting import ComboRun_Plots

import copy
import numpy as np
import h5py

class ComboRun(ComboRun_Filters, ComboRun_Plots):

    def __init__(self, fnames):

        self.fnames = fnames

        self.peak_xs = np.array([], dtype=int)
        self.peak_ys = np.array([], dtype=int)
        self.peak_is = np.array([], dtype=int)
        self.img_inds = np.array([], dtype=int)
        self.peak_area = np.array([],dtype=int)
        self.peak_circum = np.array([],dtype=int)
        self.peak_max = np.array([],dtype=int)
        self.peak_mean = np.array([],dtype=int)
        self.peak_median = np.array([],dtype=int)
        self.peak_min = np.array([],dtype=int)
        self.peak_eccen = np.array([],dtype=int)
        self.peak_disloc= np.array([],dtype=int)

        self.fname_inds = np.array([], dtype=int)

        for i, fname in enumerate(self.fnames):
            with h5py.File(fname, 'r')  as f:
                #number of peaks in each frame. almost all 0s, size about 2000
                run_npeaks = f['5_detect/n'][:]
                run_hits_loc = np.where(run_npeaks>0)[0] #both single and multi

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
                run_peak_circum = f['6_analyse/peak_circumference'][:].ravel()
                run_peak_max = f['6_analyse/peak_max'][:].ravel()
                run_peak_mean = f['6_analyse/peak_mean'][:].ravel()
                run_peak_median = f['6_analyse/peak_median'][:].ravel()
                run_peak_min = f['6_analyse/peak_min'][:].ravel()
                run_peak_disloc = f['5_detect/dislocation'][:].ravel()
                run_peak_area = f['5_detect/area'][:].ravel()
                run_peak_eccen = f['6_analyse/peak_eccentricity'][:].ravel()

                im = f['/2_process/image'][0,...]
                self.im_y, self.im_x = im.shape

            #in frames where there is no peaks, the xy positions are -1.
            xygt0_loc = np.logical_and(run_peak_xs>0, run_peak_ys>0)

            #remove the non-peaks
            run_peak_is = run_peak_is[xygt0_loc]
            run_peak_xs = run_peak_xs[xygt0_loc]
            run_peak_ys = run_peak_ys[xygt0_loc]
            run_peak_circum = run_peak_circum[xygt0_loc]
            run_peak_max = run_peak_max[xygt0_loc]
            run_peak_mean = run_peak_mean[xygt0_loc]
            run_peak_median = run_peak_median[xygt0_loc]
            run_peak_min = run_peak_min[xygt0_loc]
            run_peak_disloc = run_peak_disloc[xygt0_loc]
            run_peak_area = run_peak_area[xygt0_loc]
            run_peak_eccen = run_peak_eccen[xygt0_loc]
            ### dont need to remove non peaks from run_img_inds for some reason...

            # for some reason, this filter needs to be done seperate to the one above
            igt0_loc = run_peak_is>0
            #remove the non-peaks
            run_peak_is = run_peak_is[igt0_loc]
            run_peak_xs = run_peak_xs[igt0_loc]
            run_peak_ys = run_peak_ys[igt0_loc]
            run_peak_circum = run_peak_circum[igt0_loc]
            run_peak_max = run_peak_max[igt0_loc]
            run_peak_mean = run_peak_mean[igt0_loc]
            run_peak_median = run_peak_median[igt0_loc]
            run_peak_min = run_peak_min[igt0_loc]
            run_peak_disloc = run_peak_disloc[igt0_loc]
            run_peak_area = run_peak_area[igt0_loc]
            run_peak_eccen = run_peak_eccen[igt0_loc]
            run_img_inds = run_img_inds[igt0_loc]
            #indices that match to self.fname_in
            run_fname_inds = np.zeros(run_peak_is.size, dtype=int)+i
            self.fname_inds = np.concatenate((run_fname_inds, self.fname_inds))
            self.img_inds = np.concatenate((run_img_inds, self.img_inds))
            self.peak_xs = np.concatenate((run_peak_xs, self.peak_xs))
            self.peak_ys = np.concatenate((run_peak_ys, self.peak_ys))
            self.peak_is = np.concatenate((run_peak_is, self.peak_is))
            self.peak_circum = np.concatenate((run_peak_circum, self.peak_circum))
            self.peak_max = np.concatenate((run_peak_max, self.peak_max))
            self.peak_mean = np.concatenate((run_peak_mean, self.peak_mean))
            self.peak_median = np.concatenate((run_peak_median, self.peak_median))
            self.peak_min = np.concatenate((run_peak_min, self.peak_min))
            self.peak_disloc = np.concatenate((run_peak_disloc, self.peak_disloc))
            self.peak_area = np.concatenate((run_peak_area, self.peak_area))
            self.peak_eccen = np.concatenate((run_peak_eccen, self.peak_eccen))


        #initialize running filter
        self.filter = np.ones(len(self.peak_is)).astype(bool)


    def copy(self):
        return copy.deepcopy(self)

    def get_thumbnails(self, r=30, n_thumbs = None):

        if n_thumbs is None:
            n_thumbs = self.peak_xs.size
            
        thumbnails = np.zeros( (n_thumbs, 2*r, 2*r), dtype=np.float32)
        print('Generating thumbnails...')
        for i_peak, (peak_x, peak_y) in enumerate(zip(self.peak_xs[:n_thumbs], self.peak_ys[:n_thumbs])):
            print(f'{i_peak}\t/{n_thumbs}', end='\r')
            
            peak_x_int, peak_y_int = int(round(peak_x)), int(round(peak_y))

            if peak_x_int<r or peak_x_int>(self.im_x -r) or peak_y_int<r or peak_y_int>(self.im_y -r):
                continue
            
            
            
            fname = self.fnames[self.fname_inds[i_peak]]
            file_i = self.img_inds[i_peak]
            
            with h5py.File(fname, 'r')  as file:
                thumbnail = file['/2_process/image'][file_i, peak_y_int-r: peak_y_int+r, peak_x_int-r:peak_x_int+r]
                thumbnails[i_peak] = thumbnail.astype(np.float32)
        print(' '*50, end='\r')    
        print('Done.')
        return thumbnails
    
        
        

