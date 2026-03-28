

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class ComboRun_Plots:

    def _scatter_plot(self, X, Y, Xlabel='', Ylabel='', f=None):
        
        if f is None: f=self.filter
        plt.figure()
        plt.scatter(X[~f], Y[~f], color='r', alpha=0.25, label='Filter: False')
        plt.scatter(X[f],  Y[f], color='g', alpha=0.25, label='Filter: True')
        plt.xlabel(Xlabel)
        plt.ylabel(Ylabel)
        plt.legend()

    def _hist_plot(self, H, f=None, bins=50, range=None):
        if f is None: f=self.filter
        plt.figure()
        plt.hist(H[~f], range=range, bins=50, color='r', alpha=0.25, label='Filter: False')
        plt.hist(H[f], range=range, bins=50, color='g', alpha=0.25, label='Filter: True')
        plt.legend()

    def hist_i(self, f=None, bins=50, range=None):
        self._hist_plot(self.peak_is**(1/6), f=f, bins=bins, range=range)
    def hist_x(self, f=None, bins=50, range=None):
        self._hist_plot(self.peak_xs, f=f, bins=bins, range=range)
    def hist_y(self, f=None, bins=50, range=None):
        self._hist_plot(self.peak_ys, f=f, bins=bins, range=range)



    
    def scatter_xy(self, f=None):
        self._scatter_plot(self.peak_xs, self.peak_ys, f=f, Xlabel='x', Ylabel='y')
    def scatter_xi(self, f=None):
        self._scatter_plot(self.peak_xs, self.peak_is**(1/6), f=f, Xlabel='x', Ylabel='I$^{1/6}$')
    def scatter_yi(self, f=None):
        self._scatter_plot(self.peak_ys, self.peak_is**(1/6), f=f, Xlabel='y', Ylabel='I$^{1/6}$')


    def view_hit(self, i, f=None):
        if f is None: f=self.filter
        fname = self.fnames[self.fname_inds[f][i]]
        file_i = self.img_inds[f][i] 

        circle = patches.Circle(
            (self.peak_xs[f][i], self.peak_ys[f][i]),
            20,
            fill=False,
            linewidth=2,
            edgecolor='red'
            )
        with h5py.File(fname, 'r')  as file:
            im = file['/2_process/image'][file_i,...]
        plt.figure()
        plt.title(f'{fname}//{file_i}')
        plt.imshow(im**(1/6), extent=[0,im.shape[0], im.shape[1], 0])
        plt.gca().add_patch(circle)
        plt.xlabel("X")
        plt.ylabel("Y")
        return im










