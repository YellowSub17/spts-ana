

import numpy as np
import h5py
import matplotlib.pyplot as plt
import umap
import autoencoder
import config
import torch
from ipywidgets import interact
import ipywidgets as widgets
from sklearn.decomposition import PCA





parser = config.get_parser()
args = parser.parse_args()


all_thbn = autoencoder.concatenate_data(f'{config.DATA_DIR}/thumbnails.h5', key='thumbnails')
train_data =autoencoder.process_thumbnails(all_thbn)


data_tensor = torch.tensor(train_data, dtype=torch.float32)
device = autoencoder.get_device()

model = autoencoder.ParticleAE(latent_dim=args.latent_dim).to(device)
model.load_state_dict(torch.load(f"{config.DATA_DIR}/model_ld{args.latent_dim}.pth", map_location=device, weights_only=True))
# 3. Get the Latent Vectors
model.eval()
with torch.no_grad():
    _, latent = model(data_tensor.to(device))
    latent = latent.cpu().numpy()

# 4. Project to 2D
pca = PCA(n_components=args.PCA)
pca_coords = pca.fit_transform(latent)



reducer = umap.UMAP(n_neighbors=30, min_dist=0.1, random_state=42)
umap_coords = reducer.fit_transform(latent)



# save to h5
with h5py.File(f'{config.DATA_DIR}/coords_ld{args.latent_dim}.h5', 'w') as f:
    f['/pca_coords'] = pca_coords
    f['/umap_coords'] = umap_coords
    f['/latent_dim'] = args.latent_dim



