[Tutorials home](Tutorials.md)

# Fix a PDB structure

## First step: load required modules

```python
# Define the path to of the package
package_dir = '/path/to/MICLOT'

# Add the 'miclot' package directory to sys.path
if package_dir not in sys.path:
    sys.path.insert(0, package_dir)
```

Now you can import the package:

```python
import miclot.utilities as mcu
```

Using this way, all commands to interact with structure are callable using *mcu*. For example: mcu.C_bond()

Or you can import directly all command, so it wil be not necessary to use *mcu.*:

```python
from miclot.utilities import *
```

## 1. Add missing atoms and hydrogens & Check protonation state

This part will explain how to ...
Using PDB2PQR-APBS with Propka (or [OpenMM](https://openmm.org/))

It can be easly (and graphically) done using the PDB2PQR server: [https://server.poissonboltzmann.org/pdb2pqr](https://server.poissonboltzmann.org/pdb2pqr). Or using the script below:

**Code**

```python
pbd2pqr_parse('2bn3.pdb')
```

**Result**

Output files are :

- 2bn3_pdb2pqr.log
- 2bn3_pqr.pdb
- 2bn3.pqr



## 2. Correlate PDB file with MDTraj topology

**Code**

```python
chainID_2_chainName, chainName_2_chainID = mdtraj_chainID_2_chainName('5azz.pdb')
print("chainID_2_chainName:", chainID_2_chainName)
print("chainName_2_chainID:", chainName_2_chainID)
```

**Result**

```
chainID_2_chainName: {0: 'A', 1: 'B', 2: 'A', 3: 'B'}
chainName_2_chainID: {'A': [0, 2], 'B': [1, 3]}
```

The command also output 3 CSV files:

- 5azz_chainID_2_chainName.csv
- 5azz_chainName_2_chainID.csv
- 5azz_mdtraj_pdb_correspondance.csv


## 3. Minimize the structure

**Code**

```python
minimize_pdb('3v16.pdb')
```

**Result**

The command also output 2 files:

- 3v16_minimized.pdb
- minimization_out.csv