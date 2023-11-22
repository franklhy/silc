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
    package_data={'silc': ['data/receptor/*', 'data/tutorial/*']}
)
