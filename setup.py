from setuptools import setup, find_packages

setup(
    name='silc',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'parmed',
        'rdkit',
        'MDAnalysis',
        'vina>=1.2.5',
        'meeko>=0.5.0',
        'py3Dmol>=2.0.0',
        'openmm',
        'parmed',
        'dill',
        'openmmtools @ git+https://github.com/choderalab/openmmtools.git#egg=0.23.1',
    ],
    #package_data={'silc': ['data/receptor/*', 'data/tutorial/*']}
    package_data={'silc': ['data/*/**',]}
#                      'data/receptor/*/**',
#                      'data/tutorial/*/**',
#                      'data/residue/*/**',
#                      ]}
)
