import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display, clear_output
import numpy as np
import autoencoder
import config

class BoxSelector:
    def __init__(self, coords, flags):
        """
        Initializes the interactive UMAP bounding box selector.
        
        Parameters:
        -----------
        umap_coords : numpy.ndarray
            2D array of coordinates, shape (N, 2)
        flags : numpy.ndarray or list
            Boolean array or indices indicating 'red' vs 'green' points
        """
        self.coords = coords
        self.flags = flags

        # 1. Initialize sliders with default values matching your viewport
        self.x1_s = widgets.FloatSlider(value=0, min=self.coords[:,0].min(), max=self.coords[:,0].max(), step=0.1, description='x1:')
        self.x2_s = widgets.FloatSlider(value=1, min=self.coords[:,0].min(), max=self.coords[:,0].max(), step=0.1, description='x2:')
        self.y1_s = widgets.FloatSlider(value=0, min=self.coords[:,1].min(), max=self.coords[:,1].max(), step=0.1, description='y1:')
        self.y2_s = widgets.FloatSlider(value=1, min=self.coords[:,1].min(), max=self.coords[:,1].max(), step=0.1, description='y2:')

        thumbnails = autoencoder.concatenate_data(f'{config.DATA_DIR}/thumbnails.h5')
        self.proc = autoencoder.process_thumbnails(thumbnails)

        # 2. Create the Button and Text Output UI components
        self.print_btn = widgets.Button(description="Print Values", button_style='info', icon='print')
        self.show_btn = widgets.Button(description="Show random peaks", button_style='info', icon='print')
        self.text_output = widgets.Output()
        self.plot_output = widgets.Output()
        self.plot2_output = widgets.Output()

        # 3. Connect actions
        self.print_btn.on_click(self._print_box_coords)
        self.show_btn.on_click(self._show_peaks)
        for slider in [self.x1_s, self.x2_s, self.y1_s, self.y2_s]:
            slider.observe(self._update_plot, names='value')

    def _update_plot(self, change=None):
        """Internal method to update the matplotlib figure."""
        with self.plot_output:
            clear_output(wait=True)
            
            plt.figure(figsize=(6, 6))
            
            # Draw scatter points
            plt.scatter(self.coords[self.flags, 0], self.coords[self.flags, 1], c='green', alpha=0.3, s=10)
            plt.scatter(self.coords[~self.flags, 0], self.coords[~self.flags, 1], c='red', alpha=0.3, s=10)
            
            # Draw selection box using correct data coordinates
            plt.vlines(x=[self.x1_s.value, self.x2_s.value], ymin=self.y1_s.value, ymax=self.y2_s.value, colors='blue')
            plt.hlines(y=[self.y1_s.value, self.y2_s.value], xmin=self.x1_s.value, xmax=self.x2_s.value, colors='blue')

      
            
            plt.show()

    def _show_peaks(self, b):
        xmin = min(self.x1_s.value, self.x2_s.value)
        xmax = max(self.x1_s.value, self.x2_s.value)
        ymin = min(self.y1_s.value, self.y2_s.value)
        ymax = max(self.y1_s.value, self.y2_s.value)
         # Use bitwise & for clean, multi-condition filtering
        x_in_box = (self.coords[:, 0] > xmin) & (self.coords[:, 0] < xmax)
        y_in_box = (self.coords[:, 1] > ymin) & (self.coords[:, 1] < ymax)
        # Combine them to get the final indices
        inds = np.where(x_in_box & y_in_box)[0]

      
        with self.plot2_output:
            clear_output(wait=True)
            
            fig, axs = plt.subplots(4,4,figsize=(6, 6))
            for ax in axs.flatten():
                i = np.random.choice(inds)
                ax.imshow(self.proc[i, 0,:,:])

            plt.show()
            




   



    def _print_box_coords(self, b):
        """Internal method triggered by clicking the print button."""
        xmin = min(self.x1_s.value, self.x2_s.value)
        xmax = max(self.x1_s.value, self.x2_s.value)
        ymin = min(self.y1_s.value, self.y2_s.value)
        ymax = max(self.y1_s.value, self.y2_s.value)
         # Use bitwise & for clean, multi-condition filtering
        x_in_box = (self.coords[:, 0] > xmin) & (self.coords[:, 0] < xmax)
        y_in_box = (self.coords[:, 1] > ymin) & (self.coords[:, 1] < ymax)
        # Combine them to get the final indices
        inds = np.where(x_in_box & y_in_box)[0]
  
        with self.text_output:
            clear_output()
            print("--- Current Box Coordinates ---")
            print(f'xmin={xmin:.2f}, xmax={xmax:.2f}, ymin={ymin:.2f}, ymax={ymax:.2f}\n')

        # fig, axs = plt.subplots(4,4,figsize=(6,6))
        # for ax in axes.flatten():

        

            # i = np.random.randint(0, proc.shape[0])
            # ax.imshow(proc[i, 0,:,:])

    def show(self):
        """Renders the complete widget control interface into the notebook."""
        slider_box = widgets.VBox([self.x1_s, self.x2_s, self.y1_s, self.y2_s, self.print_btn, self.show_btn])
        
        # Display layout panels
        display(widgets.HBox([slider_box, self.text_output]))
        # display(widgets.HBox([self.plot_output, self.plot2_output]))
        display(self.plot_output)
        display(self.plot2_output)
        
        # Trigger initial plot baseline
        self._update_plot()
