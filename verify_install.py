import sys
import os
import torch
import numpy as np

# 1. Setup Paths (The code needs to know where the sub-folders are)
# This mimics what the Colab notebook did
sys.path.append(os.path.join(os.getcwd(), "boltz_ph"))
sys.path.append(os.path.join(os.getcwd(), "LigandMPNN"))

print("--- DIAGNOSTICS ---")
print(f"Python Version: {sys.version.split()[0]}")
print(f"Numpy Version: {np.__version__}")
print(f"PyTorch Version: {torch.__version__}")

# 2. Check GPU
print("\n--- GPU CHECK ---")
if torch.cuda.is_available():
    print(f"✅ CUDA is available!")
    print(f"   Device: {torch.cuda.get_device_name(0)}")
    print(f"   VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.2f} GB")
else:
    print("❌ CUDA NOT DETECTED. The model will run very slowly on CPU.")

# 3. Check Imports
print("\n--- LIBRARY IMPORT CHECK ---")
try:
    from boltz.main import BoltzInference
    print("✅ Boltz imported successfully")
except ImportError as e:
    print(f"❌ Failed to import Boltz: {e}")

try:
    import model_utils
    print("✅ LigandMPNN imported successfully")
except ImportError as e:
    print(f"❌ Failed to import LigandMPNN: {e}")

print("\n-------------------")
if torch.cuda.is_available():
    print("READY TO HUNT! 🏹")
else:
    print("Completed with warnings (No GPU).")

