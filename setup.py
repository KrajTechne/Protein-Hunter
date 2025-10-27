from setuptools import setup, find_packages
from setuptools.command.install import install
from pathlib import Path
import subprocess
import sys
import os

class PostInstallCommand(install):
    """Custom post-installation for downloading Boltz weights."""
    def run(self):
        install.run(self)
        
        print("🩹 Fixing NumPy/Numba compatibility...")
        subprocess.check_call([
            sys.executable, '-m', 'pip', 'install', 
            '--upgrade', 'numpy>=1.24,<1.27', 'numba'
        ])
        
        print("⬇️  Downloading Boltz weights and dependencies...")
        try:
            from boltz.main import download_boltz2
            cache = Path('~/.boltz').expanduser()
            cache.mkdir(parents=True, exist_ok=True)
            download_boltz2(cache)
            print('✅ Boltz weights downloaded successfully!')
        except Exception as e:
            print(f"⚠️  Warning: Could not download Boltz weights: {e}")
            print("   You may need to run this manually after installation.")
        
        # Make DAlphaBall executable
        try:
            dalphaball = Path(__file__).parent / "boltz" / "utils" / "DAlphaBall.gcc"
            if dalphaball.exists():
                os.chmod(dalphaball, 0o755)
                print("✅ DAlphaBall.gcc set as executable")
        except Exception as e:
            print(f"⚠️  Warning: Could not chmod DAlphaBall.gcc: {e}")

setup(
    name='protein-hunter',
    version='0.1.0',
    description='Protein Hunter',
    long_description=open('README.md').read() if Path('README.md').exists() else 'Protein Hunter',
    long_description_content_type='text/markdown',
    author='Yehlin Cho',
    author_email='yehlin@mit.edu',
    url='https://github.com/yehlincho/Protein-Hunter',
    packages=find_packages(include=['boltz*']),
    include_package_data=True,
    python_requires='>=3.10,<3.11',
    install_requires=[
        # Core dependencies
        'matplotlib',
        'seaborn',
        'prody',
        'tqdm',
        'PyYAML',
        'requests',
        'pypdb',
        'py3Dmol',
        'logmd==0.1.45',
        'ml_collections',
        # NumPy/Numba will be upgraded in post-install
        'numpy>=1.23',
        'numba',
        # Jupyter support
        'ipykernel',
    ],
    extras_require={
        'pyrosetta': [
            'pyrosettacolabsetup',
            'pyrosetta-installer',
        ],
        'dev': [
            'pytest',
            'pytest-cov',
            'black',
            'flake8',
        ],
    },
    cmdclass={
        'install': PostInstallCommand,
    },
    entry_points={
        'console_scripts': [
            'boltz-setup-jupyter=boltz_design.scripts:setup_jupyter_kernel',
        ],
    },
    package_data={
        'boltz': ['utils/DAlphaBall.gcc'],
    },
    zip_safe=False,
)