[Tutorials home](Tutorials.md)

# Tutorial: How to identify non-bonded interaction ?

This tutorial give you working example on how to identify non-bonded interaction between two residue, or indide one residue: only for C5 H-bond.

## First step: load required modules

```python
import ***
```


## 0. C5 hydrogen bond

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

### 2.1. No interaction

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

### 2.2. Hydrophobic interaction

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

### 2.3. Hydrophobe/Hydrophile clash

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

### 2.4. Hydrophobe/Hydrophile without clash

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

### 3.1. Charge clash

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

### 3.2. Charge repulsion

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


## 5. Hydrogen bond

> [!IMPORTANT]
> Remember, as demonstrated in the previous section, that the chosen method influences the output. 

### 5.1. Simple usage example

**Code**

```python
pdb_file = 'Hbond_1.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = hydrogen_bond(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_angle)
print(interaction.get_distance)
print(interaction.get_subtype)
```

**Result**

```
True
[[12 18  6]
 [17 22  6]]
[[151.36641, array([12, 18,  6])], [156.7908, array([17, 22,  6])]]
[[3.152380883693695, array([12, 18,  6])], [3.1175026297569275, array([17, 22,  6])]]
[['regular', array([12, 18,  6])], ['regular', array([17, 22,  6])]]
```


### 5.2. Effect of changing parameters of the "baker_hubbard"

**Code**

```python
pdb_file = 'Hbond_2.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = hydrogen_bond(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_angle)
print(interaction.get_distance)
print(interaction.get_subtype)
```

**Result**

```
False
[]
[]
[]
[]
```

Now, if we decrese the angle cutoff from 150° (Default) to 120°, the command detect an H-bond.

**Code**

```python
pdb_file = 'Hbond_2.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = hydrogen_bond(traj, 0,1, angle_cutoff=120)
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_angle)
print(interaction.get_distance)
print(interaction.get_subtype)
```

**Result**

```
True
[[ 7 14 21]
 [ 7 14 28]]
[[122.01538, array([ 7, 14, 21])], [145.35005, array([ 7, 14, 28])]]
[[3.183967173099518, array([ 7, 14, 21])], [3.0763807892799377, array([ 7, 14, 28])]]
[['regular', array([ 7, 14, 21])], ['regular', array([ 7, 14, 28])]]
```




## 6. van der Waals

### 6.1. Simple usage

**Code**

```python
pdb_file = 'vdw.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = van_der_waals(traj, 1,2)
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_number_contacts)
print(interaction.get_interface)
```

**Result**

```
True
([2.7280741930007935, 2.323926091194153, 2.363763749599457, 3.7214794754981995, 2.4790769815444946, 2.4385419487953186, 3.139965832233429, 2.0333456993103027, 2.530612051486969, 3.2540124654769897, 2.2418467700481415, 2.785297930240631, 3.1580671668052673, 2.5184500217437744, 3.466078042984009, 3.0537447333335876], [[ARG137-N, GLY138-N], [ARG137-N, GLY138-H], [ARG137-CA, GLY138-N], [ARG137-CA, GLY138-CA], [ARG137-CA, GLY138-H], [ARG137-C, GLY138-CA], [ARG137-C, GLY138-C], [ARG137-C, GLY138-H], [ARG137-C, GLY138-HA3], [ARG137-C, GLY138-HA2], [ARG137-O, GLY138-N], [ARG137-O, GLY138-CA], [ARG137-O, GLY138-C], [ARG137-O, GLY138-HA3], [ARG137-CB, GLY138-N], [ARG137-HA, GLY138-N]])
16
24.29358959197998
```



### 6.2. Effect of inoring hydrogens

Ignoring hydrogens during this analysis decrease the number of contact and the interface between the two residues.

**Code**

```python
pdb_file = 'vdw.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = van_der_waals(traj, 1,2 , set_hydrogen=False)
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_number_contacts)
print(interaction.get_interface)
```

**Result**

```
True
([2.7280741930007935, 2.363763749599457, 3.7214794754981995, 2.4385419487953186, 3.139965832233429, 2.2418467700481415, 2.785297930240631, 3.1580671668052673, 3.466078042984009], [[ARG137-N, GLY138-N], [ARG137-CA, GLY138-N], [ARG137-CA, GLY138-CA], [ARG137-C, GLY138-CA], [ARG137-C, GLY138-C], [ARG137-O, GLY138-N], [ARG137-O, GLY138-CA], [ARG137-O, GLY138-C], [ARG137-CB, GLY138-N]])
9
22.14486598968506
```



### 6.3. Effect of reducing the `distance_toerance` parameter

Playing with the distance tolerance *N* can increase or decrease the number of contact and the interface between the two residues.

**Code**

```python
pdb_file = 'vdw.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = van_der_waals(traj, 1,2 , distance_tolerance=0)
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_number_contacts)
print(interaction.get_interface)
```

**Result**

```
True
([2.7280741930007935, 2.323926091194153, 2.363763749599457, 2.4790769815444946, 2.4385419487953186, 3.139965832233429, 2.0333456993103027, 2.530612051486969, 2.2418467700481415, 2.785297930240631, 3.1580671668052673, 2.5184500217437744], [[ARG137-N, GLY138-N], [ARG137-N, GLY138-H], [ARG137-CA, GLY138-N], [ARG137-CA, GLY138-H], [ARG137-C, GLY138-CA], [ARG137-C, GLY138-C], [ARG137-C, GLY138-H], [ARG137-C, GLY138-HA3], [ARG137-O, GLY138-N], [ARG137-O, GLY138-CA], [ARG137-O, GLY138-C], [ARG137-O, GLY138-HA3]])
12
24.29358959197998
```



### 6.4. No interaction *between residues indices 0 and 1*

**Code**

```python
pdb_file = 'vdw.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = van_der_waals(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_number_contacts)
print(interaction.get_interface)
```

**Result**

```
False
([], [])
0
None
```




## 7. Amino-$\pi$

### 7.1. Simple usage

**Code**

```python
pdb_file = '7veg.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = amino_pi(traj, 12,54) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
True
3.5438241270201325
109.76115407329473

```


### 7.2. No interaction

**Code**

```python
pdb_file = '7veg.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = amino_pi(traj, 34,10) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
False
7.389154148287681
62.269396599020865
```





<!--- TEMPLATE
## 4. Salt bridge

### 3.1. Simple usage

**Code**

```python
```

**Result**

```
```
-->
