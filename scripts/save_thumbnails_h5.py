import numpy as np
import matplotlib.pyplot as plt
import comborun
from config import generate_full_paths, ps_data_ranges
import h5py


file_path = '../data/thumbnails.h5'


with h5py.File(file_path, "w") as f:
    for group in ps_data_ranges.keys():
        print(group)
        cbr1 = comborun.ComboRun(generate_full_paths(ps_data_ranges[group], 25))
        cbr2 = comborun.ComboRun(generate_full_paths(ps_data_ranges[group], 7))
        cbr1.filter_focused(cbr2)
        thumbnails = cbr1.get_thumbnails()
        # --- H5PY WRITE START ---

        grp = f.create_group(group)
        grp.create_dataset("thumbnails", data=thumbnails, compression="gzip", chunks=(1, 60, 60))
        grp.create_dataset("flags", data=cbr1.filter[:] )
        grp.create_dataset("xs", data=cbr1.peak_xs[:] )
        grp.create_dataset("ys", data=cbr1.peak_ys[:])
        grp.create_dataset("is", data=cbr1.peak_is[:])
        grp.create_dataset("circum", data=cbr1.peak_circum[:])
        grp.create_dataset("max", data=cbr1.peak_max[:])
        grp.create_dataset("mean", data=cbr1.peak_mean[:])
        grp.create_dataset("median", data=cbr1.peak_median[:])
        grp.create_dataset("min", data=cbr1.peak_min[:])
        grp.create_dataset("disloc", data=cbr1.peak_disloc[:])
        grp.create_dataset("area", data=cbr1.peak_area[:])
        grp.create_dataset("eccen", data=cbr1.peak_eccen[:])

        
        # --- H5PY WRITE END ---

print(f"All groups successfully written to {file_path}")



# for group in ps_data_ranges.keys():
    # print(group)

    # cbr1 = comborun.ComboRun(generate_full_paths(ps_data_ranges[group],25))
    # cbr2 = comborun.ComboRun(generate_full_paths(ps_data_ranges[group],7))
    # cbr1.filter_focused(cbr2)

    # thumbnails = cbr1.get_thumbnails()

    # flags = cbr1.filter[:]


    # output[group] = {
            # 'thumbnails': thumbnails,
            # 'flags': flags
            # }







