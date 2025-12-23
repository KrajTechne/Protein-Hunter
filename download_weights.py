import sys
import os
from pathlib import Path

# 1. Add the specific source folder to path (just like Colab)
# This is critical because the code is structured in a non-standard way
sys.path.insert(0, os.path.join(os.getcwd(), "boltz_ph"))

try:
    from boltz.main import download_boltz2
    print("✅ Found downloader function.")
except ImportError:
    # Fallback: sometimes the import path varies slightly
    sys.path.insert(0, os.path.join(os.getcwd(), "boltz_ph", "src"))
    from boltz.main import download_boltz2
    print("✅ Found downloader function (via src).")

# 2. Set the destination to the default location Boltz expects (~/.boltz)
# If we download to ".", the main script might not find it later.
# The error you saw looked for: ~/.boltz/boltz2_conf.ckpt
home_path = Path.home() / ".boltz"
os.makedirs(home_path, exist_ok=True)

print(f"⬇️ Downloading weights and CCD to: {home_path}")
download_boltz2(home_path)

print("🎉 Download Complete! You can now run the design script.")
