
import numpy as np
import matplotlib.pyplot as plt
import h5py

import config
import comborun
import autoencoder
from config import generate_full_paths, ps_data_ranges, get_parser

from torch.utils.data import DataLoader, TensorDataset
import torch
import torch.nn as nn



parser = get_parser()
args = parser.parse_args()


all_thbn = autoencoder.concatenate_data(f'{config.DATA_DIR}/thumbnails.h5')

train_data = autoencoder.process_thumbnails(all_thbn)







# Convert to PyTorch Tensors
# Filter out the empty (zero) thumbnails first!
valid_mask = train_data.sum(axis=(1, 2, 3)) > 0
data_tensor = torch.tensor(train_data[valid_mask], dtype=torch.float32)

dataset = TensorDataset(data_tensor)
loader = DataLoader(dataset, batch_size=32, shuffle=True)

# Initialize

device = autoencoder.get_device()

model = autoencoder.ParticleAE(latent_dim=args.latent_dim).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
criterion = nn.MSELoss()



# Loop
losses = []
for epoch in range(args.epochs):
    total_loss = 0
    for batch in loader:
        img = batch[0].to(device)
        # Forward
        recon, _ = model(img)
        loss = criterion(recon, img)
        # Backward
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    losses.append(round(total_loss/len(loader), 6))
    if epoch % 10 == 0:
        print(f"Epoch {epoch}, Loss: {total_loss/len(loader):.6f}")


model_path = f"{config.DATA_DIR}/model_ld{args.latent_dim}.pth"

# Save the weights
torch.save(model.state_dict(), model_path)
print(f"Model saved to {model_path}")


# plt.figure()
# plt.plot(losses)
# plt.show()
