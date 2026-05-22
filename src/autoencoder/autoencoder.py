import torch
import torch.nn as nn

import numpy as np

from torch.utils.data import DataLoader, TensorDataset



class ParticleAE(nn.Module):
    def __init__(self, latent_dim=32):
        super(ParticleAE, self).__init__()
        
        # Encoder: 1x60x60 -> latent_dim
        self.encoder = nn.Sequential(
            nn.Conv2d(1, 16, 3, stride=2, padding=1),  # -> 16x30x30
            nn.ReLU(),
            nn.Conv2d(16, 32, 3, stride=2, padding=1), # -> 32x15x15
            nn.ReLU(),
            nn.Flatten(),
            nn.Linear(32 * 15 * 15, latent_dim)
        )
        
        # Decoder: latent_dim -> 1x60x60
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 32 * 15 * 15),
            nn.Unflatten(1, (32, 15, 15)),
            nn.ConvTranspose2d(32, 16, 3, stride=2, padding=1, output_padding=1), # -> 16x30x30
            nn.ReLU(),
            nn.ConvTranspose2d(16, 1, 3, stride=2, padding=1, output_padding=1),  # -> 1x60x60
            nn.Sigmoid() # Keeps output pixels between 0 and 1
        )

    def forward(self, x):
        latent = self.encoder(x)
        reconstruction = self.decoder(latent)
        return reconstruction, latent


def process_thumbnails(thumbnails, sigma=15):
    # 1. Create the Gaussian Mask (Same size as a single thumbnail)
    y, x = np.ogrid[-30:30, -30:30]
    mask = np.exp(-(x**2 + y**2) / (2 * sigma**2))
    
    # 2. Background Subtraction (Per-thumbnail median)
    # We calculate the median for each 60x60 block and subtract it
    medians = np.median(thumbnails, axis=(1, 2), keepdims=True)
    processed = thumbnails - medians
    processed[processed < 0] = 0 # Clip negative noise
    
    # 3. Local Normalization
    # Find max of each thumbnail
    max_vals = processed.max(axis=(1, 2), keepdims=True)
    
    # Use np.where to avoid division by zero on empty frames
    # If max is 0, we just keep it as 0.
    processed = np.divide(processed, max_vals, out=np.zeros_like(processed), where=max_vals!=0)
    
    # 4. Apply Gaussian Mask to the whole stack
    # NumPy will automatically broadcast the (60, 60) mask across the (2000, 60, 60) array
    processed *= mask
    processed = processed[:, np.newaxis, :, :]
    return processed.astype(np.float32)