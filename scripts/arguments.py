

import argparse



def get_parser():
    parser = argparse.ArgumentParser(description="Shared Autoencoder Configurations")

    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--latent-dim", type=int, default=32)
    parser.add_argument("--PCA", type=int, default=2)

    return parser
