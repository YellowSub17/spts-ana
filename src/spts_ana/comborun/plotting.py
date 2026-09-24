

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import h5py

class ComboRun_Plots:

    def _scatter_plot(self, X, Y, Xlabel='', Ylabel='', f=None, fig=None, axes=None):
        
        if f is None: f=self.filter
        if fig is None:
            fig, axes = plt.subplots(1,1)
        axes.scatter(X[~f], Y[~f], color='r', alpha=0.25, label='Filter: False')
        axes.scatter(X[f],  Y[f], color='g', alpha=0.25, label='Filter: True')
        axes.set_xlabel(Xlabel)
        axes.set_ylabel(Ylabel)
        axes.legend()
        

    def _hist_plot(self, H, Xlabel='', Ylabel='', f=None, bins=50, range=None, fig=None, axes=None):
        if f is None: f=self.filter
        if fig is None:
            fig, axes = plt.subplots(1,1)
        axes.hist(H[~f], range=range, bins=bins, color='r', alpha=0.25, label='Filter: False')
        axes.hist(H[f], range=range, bins=bins, color='g', alpha=0.25, label='Filter: True')
        axes.set_xlabel(Xlabel)
        axes.set_ylabel(Ylabel)
        axes.legend()   

    def hist_i(self, f=None, bins=50, range=None, fig=None, axes=None):
        self._hist_plot(self.peak_is**(1/6), Xlabel='I$^{1/6}$', Ylabel='Freq.', f=f, bins=bins, range=range,fig=fig,axes=axes)
        
    def hist_x(self, f=None, bins=50, range=None, fig=None, axes=None):
        self._hist_plot(self.peak_xs, Xlabel='X', Ylabel='Freq.', f=f, bins=bins, range=range,fig=fig,axes=axes)
        
    def hist_y(self, f=None, bins=50, range=None, fig=None, axes=None):
        self._hist_plot(self.peak_ys, Xlabel='Y', Ylabel='Freq.', f=f, bins=bins, range=range,fig=fig,axes=axes)



    
    def scatter_xy(self, f=None, fig=None, axes=None):
        self._scatter_plot(self.peak_xs, self.peak_ys, f=f, Xlabel='x', Ylabel='y', fig=fig,axes=axes)
        plt.axis('equal')
        
    def scatter_xi(self, f=None, fig=None, axes=None):
        self._scatter_plot(self.peak_xs, self.peak_is**(1/6), f=f, Xlabel='X', Ylabel='I$^{1/6}$',fig=fig,axes=axes)
        
    def scatter_yi(self, f=None, fig=None, axes=None):
        self._scatter_plot(self.peak_ys, self.peak_is**(1/6), f=f, Xlabel='Y', Ylabel='I$^{1/6}$',fig=fig,axes=axes)


    def view_hit(self, i, f=None, fig=None, axes=None):
        if fig is None:
            fig, axes = plt.subplots(1,1)
        if f is None: f=self.filter
        fname = self.fnames[self.fname_inds[f][i]]
        file_i = self.img_inds[f][i] 

        circle = patches.Circle(
            (self.peak_xs[f][i], self.peak_ys[f][i]),
            #(100, 200),
            20,
            fill=False,
            linewidth=1,
            edgecolor='red'
            )
        with h5py.File(fname, 'r')  as file:
            im = file['/2_process/image'][file_i,...]
        plt.title(f'{fname}//{file_i}')
        axes.imshow(im**(1/6))#, extent=[0,im.shape[0], im.shape[1], 0])
        axes.add_patch(circle)
        axes.set_xlabel("X")
        axes.set_ylabel("Y")
        return im










