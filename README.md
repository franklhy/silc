# Silc: <ins>S</ins>creening <ins>I</ins>nter<ins>L</ins>ocked <ins>C</ins>hemistry

## Dependencies

Silc requires a number of Python packages. In order to install and run silc smoothly, the following requirements and dependencies should be fulfilled:

- python3 (tested for 3.8.18)
- RDKit: https://github.com/rdkit/rdkit
- vina: Python binder of Autodock Vina, see https://autodock-vina.readthedocs.io/en/latest/installation.html#python-bindings
- Meeko: https://github.com/forlilab/Meeko
- py3Dmol (optinoal, but it is used in tutorial): https://pypi.org/project/py3Dmol/
- AmberTool23: https://ambermd.org/GetAmber.php#ambertools
- openmm: http://docs.openmm.org/
- parmed: https://parmed.github.io/ParmEd/html/index.html
- MDAnalysis: https://www.mdanalysis.org/pages/installation_quick_start/

## Install this package

Run the following in the package root directory
```
pip install [-e] .
```
where the optional `-e` flag specifies that we want to install in *editable* mode, which means that when we edit the files in our package we do not need to re-install the package before the changes come into effect. You will need to either restart Python or reload the package though!

## Install all dependent packages

#### on midway2
Run the following:
```
mamba create -n pysages -c conda-forge python=3.9
mamba activate pysages
mamba install -c conda-forge cuda-version=11.2 cudatoolkit-dev=11.2 cudnn=8.2 cutensor nccl cupy  # CUDA deps
mamba install -c conda-forge cython numba  # compiler tools
mamba install -c conda-forge typing-extensions plum-dispatch  # typing tools
mamba install -c conda-forge importlib-metadata  # python extras

# Install jax
mamba install -c conda-forge ml_dtypes opt-einsum scipy zipp  # jax deps
pip install jaxlib==0.4.2+cuda11.cudnn82 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install jax==0.4.2 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Install openmm-dlext, installs openmm as well
mamba install -c conda-forge cuda-version=11.2 openmm-dlext

# Install pysages
pip install git+https://github.com/gustavor101/PySAGES.git@funnel_abf
```

## Tutorial

In the `tutorial` folder, I demonstrate the usage of this package. You can simply run tutorial_*.ipynb in sequence using Jupyter Notebook.
