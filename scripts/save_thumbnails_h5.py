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
        flags = cbr1.filter[:]
        # --- H5PY WRITE START ---
        # 1. Create a group in the H5 file matching your experimental group name
        grp = f.create_group(group)
        # 2. Save the 3D thumbnail array with gzip compression to save space
        grp.create_dataset(
            "thumbnails",
            data=thumbnails,
            compression="gzip",
            chunks=(1, 60, 60)  # Optimizes reading single images later
        )
        # 3. Save the corresponding 1D boolean array
        grp.create_dataset(
            "flags",
            data=flags
        )
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







