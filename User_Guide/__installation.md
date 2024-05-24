[User Guide home](Manual.md)
# Installation guide

This script is written in [Python](https://www.python.org/), with the **version 3.10.13** packaged by conda-forge.

## 1. Install Python

We recommend using a Python environment with [Conda](https://docs.conda.io/projects/conda/en/stable/index.html). You can install it with [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) or [Anaconda](https://www.anaconda.com/download/).

## 2. Install dependencies

| Name   | Version | Usage | Description | Link |
| ------ | ------- | ----- | ----------- | ---- |
| Numpy  |||||
| MDTraj |||||
| Pandas |||||
| scikit-spatial |||||
| OpenMM |||||

### 2.1 Simple way: use the *miclot.yml* file

It use an environment *.yml* file to install all the dependencies.
More technical informations concerning this file type in Conda is available [here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

```bash
conda env create -f miclot.yml
```

### 2.2 Complex way

You can do it manually, but be sure to use the same version of each dependancy and python.
For Conda, you can find information on how to install a specific version of a package on [Anaconda website](https://anaconda.org/).


## 3. Download from GitHub

### 3.1. Clone the repository

```bash
git clone https://github.com/miclot/miclot.git
```

or with the GitHub CLI

```bash
gh repo clone TMiclot/MICLOT
```

### 3.2. Download the repository as *.zip* file

On the main page, clik on the `<> Code` buttun, then on `Download ZIP`.


## 4. Using *MICLOT* in a Python script

```python
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
# Define the path to of the package
package_dir = '/home/user/analysis/MICLOT/'

# Add the 'miclot' package directory to sys.path
if package_dir not in sys.path:
    sys.path.insert(0, package_dir

# import function and classes from MICLOT
import miclot.interactions as mci
from miclot.coulomb_lj import coulomb_lj

# Load the PDB file in MDTraj
pdb_file = 'MyStructure.pdb'
trajectory = md.load(pdb_file, top=pdb_file)

# Remove solvent (water, ions, ...)
trajectory.remove_solvent(inplace=True)

# Calculate Coulomb and LJ energies for residues indeces 0 and 3
energy = coulomb_lj(trajectory, 0,3 ,'/home/user/analysis/ForceField/protein.ff14SB.xml', frame=frame)

print(energy.get_energy[0]) #isolate total energy C + LJ
print(energy.get_energy_coulomb[0]) #isolate total C
print(energy.get_energy_LJ[0]) #isolate total LJ

# Check if the residues 0 and 3 perform aromatic-aromatic interaction
interaction = mci.aromatic_aromatic(trajectory, 0, 3)

print(interaction.check_interaction)
```

