from setuptools import setup, find_packages

setup(
    name='silc',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'rdkit',
        'vina>=1.2.5',
        'meeko>=0.5.0',
        'py3Dmol>=2.0.0',
        'openmm'
    ],
    #package_data={'silc': ['data/receptor/*', 'data/tutorial/*']}
    package_data={'silc': [
                      'data/receptor/*',
                      'data/tutorial/amber_charge/*',
                      'data/tutorial/build_molecules/*',
                      'data/tutorial/build_molecules/BRD/*',
                      'data/tutorial/build_molecules/BRD/amber_charge/*',
                      'data/tutorial/build_molecules/COR/*',
                      'data/tutorial/build_molecules/COR/amber_charge/*',
                      'data/tutorial/build_molecules/TLA_TLB/*',
                      'data/tutorial/build_molecules/TLA_TLB/amber_charge/*',
                      'data/tutorial/docking/*',
                      ]}
)
