# Silc: <ins>S</ins>creening <ins>I</ins>nter<ins>L</ins>ocked <ins>C</ins>hemistry

## Dependencies

Silc requires a number of python packages. In order to install and run silc smoothly, the following requirement and dependencies should be fulfilled:

- python3 (tested for 3.8.18)
- RDKit: https://github.com/rdkit/rdkit
- vina: Python binder of Autodock Vina, see https://autodock-vina.readthedocs.io/en/latest/installation.html#python-bindings
- Meeko: https://github.com/forlilab/Meeko
- py3Dmol: https://pypi.org/project/py3Dmol/

## Install

Run the following in the package root directory
```
pip install [-e] .
```
where the optional `-e` flag specifies that we want to install in *editable* mode, which means that when we edit the files in our package we do not need to re-install the package before the changes come into effect. You will need to either restart python or reload the package though!

## Tutorial

In the `tutorial` folder, I demonstrate the usage of this package. You can simple run tutorial.ipynb interactively using Jupyter Notebook.
