import numpy as np


class ComboRun_Filters:


    def filter_y(self,ymin=-1, ymax=1e6, update=True):
        #filter the y positiions to be within range ymin/ymax
        f = np.array(list(map(lambda ypos : ymin <= ypos <= ymax, self.peak_ys)))
        if update: self.filter = np.logical_and(f, self.filter) #update running filter
        return f
    
    def filter_x(self,xmin=-1, xmax=1e6, update=True):
         #filter the x positiions to be within range xmin/xmax
        f = np.array(list(map(lambda xpos : xmin <= xpos <= xmax, self.peak_xs)))
        if update: self.filter = np.logical_and(f, self.filter) #update running filter
        return f

    def filter_i(self, imin=-1, imax=1e6, sixth_root=True, update=True):
        sf = (1/6) if sixth_root else 1
        f = np.array(list(map(lambda inten : imin <= inten**(sf) <= imax, self.peak_is)))
        if update: self.filter = np.logical_and(f, self.filter) #update running filter
        return f
        
    
    def filter_focused(self, r2, focus_threshold=0.9, update=True):
        # r1 is the d15 window analysis r2 is the d5 window analysis.
        # r1 will have a higher intensity, because the window the peak is summed over is larger.
        # but, if the r2 window is close to the r1 intensity (within the focus_threshold), 
        # then the peak is tightly confined within the window
        # this indicates the particle is in focus
        # so, we want the particles where r1 > r2 and r2 > r1*focus_threshold

        # f must be in self's (r1's) peak order, so look up each r1 peak's matching r2 peak
        # by (run, x, y) rather than relying on the two analyses listing peaks in the same
        # order -- spts does not guarantee that, and when the orders differ the flags end up
        # attached to the wrong particles.
        # key for the dictionary: which run (fname index, same file list for r1 and r2) and xy
        xyb_dict = {}
        for fname_ind, peak_b_x, peak_b_y, peak_b_i in zip(r2.fname_inds, r2.peak_xs, r2.peak_ys, r2.peak_is):
            xyb_dict[(fname_ind, peak_b_x, peak_b_y)] = peak_b_i
        if len(xyb_dict) != len(r2.peak_is):
            print(f'Warning: {len(r2.peak_is) - len(xyb_dict)} r2 peaks share a (run, x, y) key with another peak')

        # for each "a" peak, find the "b" peak at the same position in the same run.
        # peaks with no match are not flagged as focused.
        f = np.zeros(len(self.peak_is), dtype=bool)
        n_unmatched = 0
        for n, (fname_ind, peak_a_x, peak_a_y, peak_a_i) in enumerate(zip(self.fname_inds, self.peak_xs, self.peak_ys, self.peak_is)):
            peak_b_i = xyb_dict.get((fname_ind, peak_a_x, peak_a_y))
            if peak_b_i is None:
                n_unmatched += 1
                continue
            # b is within a and 0.9*a
            f[n] = peak_a_i**(1/6) >= peak_b_i**(1/6) >= focus_threshold*peak_a_i**(1/6)
        if n_unmatched:
            print(f'Warning: {n_unmatched} of {len(self.peak_is)} r1 peaks have no matching r2 peak')

        if update: self.filter = np.logical_and(f, self.filter)
        return f