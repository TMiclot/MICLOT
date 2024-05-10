[Tutorials home](Tutorials.md)

# Get energetic properties

## 1. Coulomb & Lennard-Jones

**Code**

```python
# Force field to use
ff = 'openmmforcefields_ffxml/protein.ff14SB.xml'

# trajectory (or structure) to use
pdb_file = '3v16.pdb'

# Load molecular structure using MDTraj
traj = md.load(pdb_file, top=pdb_file)

# residues indeces
resID_A = 260
resID_B = 219

# compure LJ and coulomb
interaction_energy = coulomb_lj(traj, resID_A, resID_B, ff)

# print results
print(interaction_energy.get_energy)
print(interaction_energy.get_energy_LJ)
print(interaction_energy.get_energy_coulomb)
```

**Result**

```
(27.29836082439392, 14.20622389018442)
(-0.11023730074672013, -0.09407264655416284, -0.06197796301079388)
(27.408598125140642, 14.300296536738584, 1.301674243893409)
```



## 2. Contact-based method

**Code**

```python
compute_binding_energy('3bzd.pdb', ['A'], ['B'])
```

**Result**

```
contacts_apolar-polar 	              15.0
contacts_apolar-apolar 	              12.0
contacts_charged-polar 	              7.0
contacts_polar-polar 	              7.0
contacts_apolar-charged               6.0
contacts_charged-charged              4.0
TOTAL_contacts 	                      51.0
NIS_polar 	                          103.0
NIS_apolar 	                          74.0
NIS_charged 	                      74.0
TOTAL_NIS_types 	                  251.0
NIS_polar(%) 	                      41.035857
NIS_apolar(%) 	                      29.482072
NIS_charged(%) 	                      29.482072
temperature(C) 	                      25.0
temperature(K) 	                      298.15
DeltaG(kcal/mol) 	                  -9.37332
dissociation_constant_(M)_at_25(C)    1.33265e-07
PDB_file 	                          3bzd.pdb
chains_receptor 	                  A
chains_ligand 	                      B
```
