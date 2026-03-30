from pathlib import Path

# This locates the root relative to THIS config file
PROJ_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJ_ROOT / "data" 
