[Tutorials home](Tutorials.md)
# How to identify a given non-bonded interaction into an amino acid pair

## First step: load required modules
```python
import ***
```

## 1. C-bond
**Code**
```python
pdb_file = '3ntu_c-bond_269-261.pdb'
traj = md.load(pdb_file, top=pdb_file)
index_A = 255 #Lys261
index_B = 256 #Ala269

interaction = C_bond(traj, index_A, index_B)
print(interaction.check_interaction)
print("Atoms:", interaction.get_atoms)
print("Angle:", interaction.get_angle)
print("Distance:", interaction.get_distance)
print("Energy:", interaction.get_energy)
```
**Result**
```
True
Atoms: [{1953, 1954, 1946, 1943}]
Angle: [[[169.06029, 170.11108], {1953, 1954, 1946, 1943}]]
Distance: [[2.871696650981903, {1953, 1954, 1946, 1943}]]
Energy: [[-21.651590464607594, {1953, 1954, 1946, 1943}]]
```

## 2. Hydrophobic interaction & Hydrophobe/Hydrophile clash
### 2.1. Case: No interaction
**Code**
```python
pdb_file = 'hydrophobe_1_no_interaction.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = hydrophobic(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
```
**Result**
```
(False, False)
(8.198665976524353, 7.304799486171744)

```

### 2.2. Case: Hydrophobic interaction
**Code**
```python
pdb_file = 'hydrophobe_2.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = hydrophobic(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
```
**Result**
```
(True, True)
(8.934004306793213, 4.607366824181708)
```

### 2.3. Case: Hydrophobe/Hydrophile clash
**Code**
```python
pdb_file = 'hydrophobe_hydrophilic_clash.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = hydrophobic(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
```
**Result**
```
(True, False)
(7.733982801437378, 4.12462471880705)
```

### 2.4. Case: Hydrophobe/Hydrophile NO clash
**Code**
```python
pdb_file = 'hydrophobe_hydrophilic_NOclash.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = hydrophobic(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
```
**Result**
```
(False, False)
(8.182007074356079, 7.443699842878453)
```


***
# How to identify internal non-bonded interaction inside a amino acid
## C5 hydrogen bond
**Code**
```python
pdb_file = 'peptide_C5_Hbond.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = C5_hydrogen_bond(traj, 2)
print(interaction.check_interaction)
print(interaction.get_angle)
print(interaction.get_distance)
```
**Result**
```
True
(-149.50195, 132.94899)
2.4230431020259857
```
