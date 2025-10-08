[User Guide home](Manual.md)
# Installation guide

> [!NOTE]  
> The software has been written and tested on Linux (Ubuntu 22.04.4 LTS), it should be compatible with other OS but this has not been tested.

This script is written in [Python](https://www.python.org/), with the **version 3.10.12** packaged by conda-forge.

## 1. Install Python

We recommend using a Python environment with [Conda](https://docs.conda.io/projects/conda/en/stable/index.html). You can install it with [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) or [Anaconda](https://www.anaconda.com/download/).

## 2. Install dependencies

| Name           | Version |
| -------------- | ------- |
| numpy          | 1.24.3  |
| mdtraj         | 1.9.9   |
| pandas         | 1.5.3   |
| scikit-spatial | 7.2.0   |
| openmm         | 8.1.1   |
| pdb2pqr        | 3.6.1   |
| propka         | 3.5.1   |
| biopython      | 1.83    |
| freesasa       | 2.2.1   |

It can be done simply using this command:

```bash
conda create -n miclot -c conda-forge python=3.10.13 freesasa=2.2.1 biopython=1.83 propka=3.5.1 pdb2pqr=3.6.1 openmm=8.1.1 scikit-spatial=7.2.0 pandas=1.5.3 mdtraj=1.9.9 numpy=1.24.3
```



## 3. Download from GitHub

### 3.1. Clone the repository

```bash
git clone https://github.com/TMiclot/MICLOT.git
```

or with the GitHub CLI

```bash
gh repo clone TMiclot/MICLOT
```

### 3.2. Download the repository as *.zip* file

On the main page, clik on the `<> Code` buttun, then on `Download ZIP`.


## 4. Using *MICLOT* in a Python script

```python
import sys

# Define the path to of the package
package_dir = '/path/to/package/MICLOT/'

# Add the 'miclot' package directory to sys.path
if package_dir not in sys.path:
    sys.path.insert(0, package_dir)
```

Then you can import the package like any other:

```python
import miclot
```
    
### Example

```python
import sys

# Define the path to of the package
package_dir = '/home/user/analysis/MICLOT/'

# Add the 'miclot' package directory to sys.path
if package_dir not in sys.path:
    sys.path.insert(0, package_dir)

# import function and classes from MICLOT
import miclot.interactions as mci
from miclot.coulomb_lj import coulomb_lj

# Load the PDB file in MDTraj
pdb_file = 'MyStructure.pdb'
trajectory = md.load(pdb_file, top=pdb_file)

# Remove solvent (water, ions, ...)
trajectory.remove_solvent(inplace=True)

# Calculate Coulomb and LJ energies for residues indeces 0 and 3
energy = coulomb_lj(trajectory, 0,3 ,'/home/user/analysis/ForceField/protein.ff14SB.xml')

print(energy.get_energy[0]) #isolate total energy C + LJ
print(energy.get_energy_coulomb[0]) #isolate total C
print(energy.get_energy_LJ[0]) #isolate total LJ

# Check if the residues 0 and 3 perform aromatic-aromatic interaction
interaction = mci.aromatic_aromatic(trajectory, 0, 3)

print(interaction.check_interaction)
```

## 5. Setup MICLOT as environment for jupyter notebook

0. Create an environment only for JupyterLab (if not already created).

```bash
conda install -c conda-forge jupyterlab
```

1. Install `nb_conda_kernels` in your jupyter environment.

```bash
conda install -n jupyterlab -c conda-forge nb_conda_kernels
```

2. Install `ipykernel` in your miclot environment.

```bash
conda install -n miclot -c conda-forge ipykernel

```
