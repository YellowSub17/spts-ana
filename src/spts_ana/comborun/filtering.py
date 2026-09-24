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

        #look up dictionary that hashes the xy position of the peaks "a"
        # key for the dictionary
        # store the peak in the dictionary
        xya_dict = {}
        for peak_a_x, peak_a_y, peak_a_i in zip(self.peak_xs, self.peak_ys, self.peak_is):
            xya_key = (peak_a_x, peak_a_y) 
            xya_dict[xya_key] = peak_a_i    


        # for each "b" peak
        # what is the key in the dictionary that would be shared with "a"
        # if the key exists there is an "a" peak that matches this "b" peak
        # add the a an b peaks to the list of pairs.
        pairs = [] 
        for peak_b_x, peak_b_y, peak_b_i in zip(r2.peak_xs, r2.peak_ys, r2.peak_is): 
            xyb_key = (peak_b_x, peak_b_y) 
            if xyb_key in xya_dict:  
                pairs.append( (xya_dict[xyb_key], peak_b_i)) 

        # pairs is a list of pairs of matching peaks
        # filter the pairs so b is with a and 0.9*a
        f = np.array(list(map( lambda pair : pair[0]**(1/6) >= pair[1]**(1/6) >= focus_threshold*pair[0]**(1/6), pairs)))
        if update: self.filter = np.logical_and(f, self.filter)
        return f