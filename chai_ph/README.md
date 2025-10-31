# chai-protein-hunter

# new features
- support for multi-chain template
- template weights (can use PDE for example)
- weighted alignment during diffusion
- partial diffusion (start with initial coordinates)

# install
```bash
sudo apt-get install -y -qq aria2
```
```python
import os
import subprocess
import threading
import site

def setup_models():
    # Setup paths
    #base_url = "https://chaiassets.com/chai1-inference-depencencies"
    base_url = "http://files.ipd.uw.edu/pub/protein_hunter/chai1"
    dl_dir = f"{site.USER_SITE}/downloads"
    os.makedirs(f"{dl_dir}/models_v2", exist_ok=True)
    os.makedirs(f"{dl_dir}/esm", exist_ok=True)

    # Download files individually with appropriate settings
    downloads = [
        (f"{dl_dir}", "conformers_v1.apkl", "-x16 -s16"),
        (f"{dl_dir}/esm", "esm2/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt", "-x16 -s16"),
        (f"{dl_dir}/models_v2", "models_v2/trunk.pt", "-x16 -s16"),
        (f"{dl_dir}/models_v2", "models_v2/diffusion_module.pt", "-x16 -s16"),
        (f"{dl_dir}/models_v2", "models_v2/confidence_head.pt", "-x8"),
        (f"{dl_dir}/models_v2", "models_v2/feature_embedding.pt", ""),
        (f"{dl_dir}/models_v2", "models_v2/token_embedder.pt", ""),
        (f"{dl_dir}/models_v2", "models_v2/bond_loss_input_proj.pt", "")
    ]

    for target_dir, url_path, opts in downloads:
        url = f"{base_url}/{url_path}"
        subprocess.run(f"aria2c {opts} --dir={target_dir} {url}", shell=True, check=True)

# Start model downloads in background
print("Starting model downloads...")
download_thread = threading.Thread(target=setup_models)
download_thread.start()

# Install chai in parallel
print("Installing chai-lab...")
os.system("pip install --no-deps git+https://github.com/sokrypton/chai-lab.git \
'gemmi~=0.6.3' 'jaxtyping>=0.2.25' 'pandera>=0.24' 'antipickle==0.2.0' \
'rdkit~=2024.9.5' 'modelcif>=1.0' 'biopython>=1.83' typing_inspect \
ihm mypy_extensions equinox wadler_lindig py3Dmol")

print("Installing LigandMPNN...")
os.system("git clone https://github.com/sokrypton/LigandMPNN.git")
os.system("mkdir model_params")
os.system("bash LigandMPNN/get_model_params.sh model_params")
os.system("pip install git+https://github.com/prody/ProDy.git")
os.system("pip install ml_collections")

# Wait for downloads to complete
print("Waiting for downloads...")
download_thread.join()
print("Setup complete!")
```
