[Tutorials home](Tutorials.md)

This tutorial give you working example on how to identify non-bonded interaction between two residue, or indide one residue: only for C5 H-bond.

## First step: load required modules

```python
import ***
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





# How to identify a given non-bonded interaction into an amino acid pair

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
(False, False, None)
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
(True, True, 'hydrophobic')
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
(True, False, 'clash')
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
(False, False, None)
(8.182007074356079, 7.443699842878453)
```



## 3. Charge clash & Charge repulsion

### 3.1. Case: Charge clash

**Code**

```python
pdb_file = 'charge_clash.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = charge_clash_repulsion(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
```

**Result**

```
(True, True, 'clash')
(4.897868633270264, 3.945750892162323)
```

### 3.2. Case: Charge repulsion

**Code**

```python
pdb_file = 'charge_repulsion.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = charge_clash_repulsion(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
```

**Result**

```
(True, False, 'repulsion')
(7.883095741271973, 6.652364134788513)
```





## 4. Salt bridge

### 4.1. Case: Interaction

**Code**

```python
pdb_file = 'salt_bridge_1.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = salt_bridge(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_hbond)
```

**Result**

```
(True, True)
(7.885891795158386, 3.931843936443329)
[[ 8 19 29]]
```

### 4.2. Case: No interaction

**Code**

```python
pdb_file = 'salt_bridge_2_NOinteraction.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = salt_bridge(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_hbond)
```

**Result**

```
(False, False)
(9.199607968330383, 5.618550777435303)
[]
```

### 4.3. Effect of presence or absence of hydrogens

**Code**

```python
pdb_file = 'salt_bridge_3.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = salt_bridge(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_hbond)
```

**Result**

```
(True, True)
(4.936801493167877, 3.6193114519119263)
[[23 36  7]
 [23 36  8]]
```

Now if we use exactly the same structure, but without hydrogens, the result is diffrent.

**Code**

```python
pdb_file = 'salt_bridge_3_noH.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = salt_bridge(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_hbond)
```

**Result**

```
(True, False)
(4.936801493167877, 3.6193114519119263)
[]
```

Previously the output was `True, True`, but wihout H the result is `True, False`.
If you read the documentation for the command [salt_bridge](../User_Guide/__usage_identify_nonbonded_interactions.md#5-salt-bridge),
this result happend because the command detect only the ionic bond and not the H-bond. 

### 4.4. Effect of the method used

The method used by MDTraj to identify H-bond are diffrent, and sometime they do not produce the same output.
So it's important to know the limitations of each method.

**Code using "baker_hubbard"**

This is the default method used by the command.

```python
pdb_file = 'salt_bridge_3.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = salt_bridge(traj, 0,1, method="baker_hubbard")
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_hbond)
```

**Result**

```
(True, True)
(4.936801493167877, 3.6193114519119263)
[[23 36  7]
 [23 36  8]]
```

**Code using "wernet_nilsson"**

```python
pdb_file = 'salt_bridge_3.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = salt_bridge(traj, 0,1, method="wernet_nilsson")
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_hbond)
```

**Result**

```
(True, True)
(4.936801493167877, 3.6193114519119263)
[[23 36  8]]
```

**Code using "kabsch_sander"**

```python
pdb_file = 'salt_bridge_3.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = salt_bridge(traj, 0,1, method="kabsch_sander")
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_hbond)
```

**Result**

```
(True, False)
(4.936801493167877, 3.6193114519119263)
[]
```


<!--- TEMPLATE
## 4. Salt bridge

### 3.1. Case: Charge clash

**Code**

```python
```

**Result**

```
```
-->
