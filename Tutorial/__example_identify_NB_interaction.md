[Tutorials home](Tutorials.md)

# Tutorial: How to identify non-bonded interaction ?

This tutorial give you working example on how to identify all non-bonded interaction into a system of interest and in a pair only. Then all individial commands will be shown.



## First step: load required modules

```python
import sys

# Define the path to of the package
package_dir = '/path/to/MICLOT'

# Add the 'miclot' package directory to sys.path
if package_dir not in sys.path:
    sys.path.insert(0, package_dir)
```

Now you can import the package:

```python
import miclot.interactions as mci
```

Using this way, all commands to analyse an interaction type are callable using *mci*. For example: mci.C_bond()

Or you can import directly all command, so it wil be not necessary to use *mci.*:

```python
from miclot.interactions import *
```



## A. Identify all interfactions in a protein, or a complex

**Code**

```python
import mdtraj as md

# trajectory (or structure) to use
pdb_file = '5azz.pdb'

# Load molecular structure using MDTraj
traj = md.load(pdb_file, top=pdb_file)

# remove solvent
traj.remove_solvent(inplace=True)

# analyse the whole system
df_interaction, df_class = mci.interaction_table_whole_system(traj)
```



## B. Identify all interfactions for a pair of residues

Insted of analyse a complete protein you can also check all possibilities for a single pair of residues.

**Code**

```python
# trajectory (or structure) to use
pdb_file = '5azz.pdb'

# Load molecular structure using MDTraj
traj = md.load(pdb_file, top=pdb_file)

# remove solvent
traj.remove_solvent(inplace=True)

# analyse between residues IDs 0 and 1
df_result, df_class = mci.identify_all_interaction_pair(traj, [0,1])

print(df_result)
print(df_class)
```

**Result**

```
   residue_1_chain  residue_1_index  residue_1_resSeq residue_1_name  \
0                0                0                 1            GLY   

   residue_2_chain  residue_2_index  residue_2_resSeq residue_2_name  \
0                0                1                 2            ILE   

   residue_1_C5Hbond  residue_2_C5Hbond  ...  salt_bridge  S_pi  \
0                NaN                NaN  ...          NaN   NaN   

   S_intermediate  S_quadrupole  Se_pi  Se_intermediate  Se_quadrupole  \
0             NaN           NaN    NaN              NaN            NaN   

   sse_hbond  sse_chalcogen  van_der_waals  
0          0              0              1  

[1 rows x 51 columns]
   residue_1_index  residue_2_index residue_1_C5Hbond residue_2_C5Hbond  \
0                0                1              None              None   

                                              c_bond amino_pi arg_involved  \
0  <miclot.interactions.C_bond object at 0x72fe94...     None         None   

  aromatic_aromatic charge_aromatic charge_clash_repulsion  \
0              None            None                   None   

                                       hydrogen-bond hydrophobic  \
0  <miclot.interactions.hydrogen_bond object at 0...        None   

                                                n_pi pi_hydrogen_bond  \
0  <miclot.interactions.n_pi object at 0x72fe9012...             None   

  salt_bridge sse_aromatic                        sse_hydrogen_chalcogen_bond  \
0        None         None  <miclot.interactions.sse_hydrogen_chalcogen_bo...   

                                       van_der_waals  
0  <miclot.interactions.van_der_waals object at 0...  
```


### B.2. Recover raw information

To recover raw information of the analysis, you can:

- Directly use the 'df_class', in the same script.
- Read the *.pkl.gz* file as Pandas dataframe, at a later inspection.

As you will see each cell of the dataframe contain the output as `<miclot.interactions>` object. This means that properties can be used as if they were the normal output of an analysis command, as you can see below:

**Code**

```python
import miclot.interactions as mci
import pandas as pd

result = pd.read_pickle("class_table_whole_system.pkl.gz")

result = pd.read_pickle("/Data/4pas/4pas_class_11_50.pkl.gz")

i = result["hydrogen_bond"][0]

print(i.check_interaction)
print(i.get_distance)
```

**Result**

```
True
[[2.8605997562408447, array([205, 206, 834])]]
```




## C. Locate interaction area

### C.1. Normal usage

**Code**

```python
pdb_file = '4pas.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = van_der_waals(traj, 52,53) # identify vdw interaction
location = locate(interaction)           # locate are for all vdw contact

l = mci.locate(class_itype)
print(l.get_interaction_type)
print('='*50)
print(l.get_atoms)
print('='*50)
print(l.get_area)
```

**Result**

```
van_der_waals
[[856, 873], [856, 879], [857, 873], [857, 874], [857, 879], [858, 874], [858, 875], [858, 877], [858, 879], [858, 880], [859, 873], [859, 874], [859, 875], [859, 876], [859, 880], [860, 873], [860, 879], [865, 873], [865, 879]]
['backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-sidechain', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-sidechain', 'sidechain-sidechain', 'backbone-sidechain', 'sidechain-sidechain']
```


### C.2. Usage with a cleaned structure information file

To get the cleaned structure information file, see [Clean structure information](../User_Guide/__usage_data_analysis.md#a2-clean-structure-information)

**Code**

```python
pdb_file = '4pas.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = van_der_waals(traj, 52,53) # identify vdw interaction
location = locate(interaction)           # locate are for all vdw contact

l = mci.locate(class_itype, directory="/home/miclot/Documents/1-analyse_ppi/5-clean_data/4pas")
print(l.get_interaction_type)
print('='*50)
print(l.get_atoms)
print('='*50)
print(l.get_area)
print('='*50)
print(l.get_code)
print('='*50)
print(l.get_table_interaction)
print('='*50)
print(l.get_table_code_name)
print('='*50)
print(l.get_table_code_simple)
print('='*50)
print(l.get_table_code_complete)
```

**Result**

```
van_der_waals
==================================================
[[856, 873], [856, 879], [857, 873], [857, 874], [857, 879], [858, 874], [858, 875], [858, 877], [858, 879], [858, 880], [859, 873], [859, 874], [859, 875], [859, 876], [859, 880], [860, 873], [860, 879], [865, 873], [865, 879]]
==================================================
['backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-sidechain', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-backbone', 'backbone-sidechain', 'sidechain-sidechain', 'backbone-sidechain', 'sidechain-sidechain']
==================================================
(['BB-BB', 'BB-BB', 'BB-BB', 'BB-BB', 'BB-BB', 'BB-BB', 'BB-BB', 'HIS-BB', 'BB-BB', 'BB-BB', 'BB-BB', 'BB-BB', 'BB-BB', 'BB-BB', 'BB-BB', 'HIS-BB', 'HIS-HIS', 'HIS-BB', 'HIS-HIS'], ['BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'HIS_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'BB_H-BB_H', 'HIS_H-BB_H', 'HIS_H-HIS_H', 'HIS_H-BB_H', 'HIS_H-HIS_H'], ['BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'HIS_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'BB_H_H-BB_H_H', 'HIS_H_H-BB_H_H', 'HIS_H_H-HIS_H_H', 'HIS_H_H-BB_H_H', 'HIS_H_H-HIS_H_H'])
==================================================
pair_index  van_der_waals_backbone-backbone  van_der_waals_backbone-sidechain  \
52_53                                    14                                 3   

pair_index  van_der_waals_sidechain-sidechain  
52_53                                       2  
==================================================
                van_der_waals
pair_code_name               
BB-BB                      14
HIS-BB                      3
HIS-HIS                     2
==================================================
                                    van_der_waals
pair_code_name_secondary_structure               
BB_H-BB_H                                      14
HIS_H-BB_H                                      3
HIS_H-HIS_H                                     2
==================================================
                    van_der_waals
pair_code_complete               
BB_H_H-BB_H_H                  14
HIS_H_H-BB_H_H                  3
HIS_H_H-HIS_H_H                 2
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

### 2.4. Hydrophobe/Hydrophile repulsion

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
(True, False, 'repulsion')
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
(True, 'charged-apolar')
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




## 7. Amino $- \pi$

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




## 8. Salt bridge

### 8.1. Case: anion- $\pi$

**Code**

```python
pdb_file = 'anion_pi.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = charge_aromatic(traj, 0, 1) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'anion-Pi')
3.6034270897476572
79.49791083897665
```

### 8.2. Case: anion-quadrupole

**Code**

```python
pdb_file = 'anion_quadrupole.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = charge_aromatic(traj, 0, 1) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'anion-quadrupole')
5.097612414323615
0.2565586234471908
```

### 8.3. Case: cation-intermediate

**Code**

```python
pdb_file = 'cation_intermediate.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = charge_aromatic(traj, 0, 1) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'cation-intermediate')
4.818607454378017
35.839345180706204
```


## 9. Aromatic-Aromatic

> [!WARNING]
> In this example, the interactions were built entirely manually. As a result, some distances are much smaller than what is actually seen. Typically, a distance for Pi-Pi stacking is around 4.5 Å up to 6 Å, according to some publications.

### 9.1. Parallel staking *with an His*

**Code**

```python
pdb_file = 'aromatic-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = aromatic_aromatic(traj, 0, 1) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'parallel')
(2.5460074915990494, 0.18071826162264856)
(17.187747281098666, 85.92963671729052, 76.88224145182211, 0.16525973727270582)
```

### 9.2. Parallel staking

**Code**

```python
pdb_file = 'aromatic-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = aromatic_aromatic(traj, 4, 5) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'parallel')
(2.595882453866925, 0.016921747589102098)
(7.1782407546847935, 89.62712644895419, 82.45980025611104, 0.026769608653339105)
```


### 9.3. Parallel offset staking

**Code**

```python
pdb_file = 'aromatic-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = aromatic_aromatic(traj, 2,3) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'offset')
(3.650135066732187, 1.6735494388664918)
(7.180590776004067, 62.710309550248354, 63.487965235398406, 0.016285738031967867)
```

### 9.4. Colanar staking

**Code**

```python
pdb_file = 'aromatic-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = aromatic_aromatic(traj, 6,7) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'coplanar')
(5.26093305171583, 5.24341615361119)
(0.8632076124065122, 4.676875951018019, 3.817067170661801, 0.0947128191129849)
```

### 9.5. Y-shaped

**Code**

```python
pdb_file = 'aromatic-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = aromatic_aromatic(traj, 8,9) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'Y-shaped')
(4.3140753842702315, 0.48750438497343557)
(86.70693939174281, 83.5115593333615, 2.532489405500428, 5.081627298469144)
```

### 9.6. T-shaped

**Code**

```python
pdb_file = 'aromatic-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = aromatic_aromatic(traj, 10,11) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'T-shaped')
(4.460427485980406, 0.3986486489188704)
(86.25358515508681, 84.87236075386747, 1.3610945987459075, 3.7199442268719167)

```




## 10. $\pi-$ hydrogen bond

### 10.1. Simple usage

**Code**

```python
pdb_file = 'Pi-HBond_5pti_clean.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = pi_hbond(traj, 0,1) 
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
True
[[25, 21]]
[[2.564013026833409, 0.25852438905004804, 3.5013544225541935]]
[[149.20084404269872]]
```

### 10.2. Simple usage 2

**Code**

```python
pdb_file = 'Pi-HBond_5pti_clean.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = pi_hbond(traj, 0,2) 
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
True
[[40, 34]]
[[2.728875808687735, 0.8294750133568008, 3.6191551379494107]]
[[143.3766091706363]]
```




## 11. $n \rightarrow \pi^*$

### 11.1. Simple usage

**Code**

```python
pdb_file = 'n_pi_bond_1k4d_res75-79.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = n_pi(traj, 1,2) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'regular')
(3.17790687084198, 4.110865294933319)
(106.11885, 99.668915)
```

### 11.2. No interaction when changing `angle_tolerance`

**Code**

```python
pdb_file = 'n_pi_bond_1k4d_res75-79.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = n_pi(traj, 0,1, angular_tolerance=5) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(False, False)
(3.051159381866455, 4.092892408370972)
(116.220634, 96.67473)
```





## 12. Chalcogen bond & S/Se mediated H-bond

### 12.1. Chalcogen bond with Se

**Code**

```python
pdb_file = 'chalcogen_SE_2gtq_withCONNECT.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = sse_hydrogen_chalcogen_bond(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'chalcogen')
([[15, 7]], [])
([3.1346604228019714], [])
([[92.7476430482533, 57.253210049674934, 122.70444577007973, 3.2657869555555123]], [])
```

### 12.2. Chalcogen bond with S

**Code**

```python
pdb_file = 'chalcogen_S_1uwc.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = sse_hydrogen_chalcogen_bond(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'chalcogen')
([[13, 5]], [])
([3.257075250148773], [])
([[87.4338907021319, 50.28551360663397, 129.6667818270496, 3.334484671896547]], [])
```

### 12.3. H-bond

**Code**

```python
pdb_file = 'hbond_S_1f9z.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = sse_hydrogen_chalcogen_bond(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'h-bond')
([], [[6, 13]])
([], [3.206118941307068])
([], [[12.765337911253035, 164.40883574237162, 77.71173636759816, 93.4844791253299]])
```

### 12.4. H-bond, changing the `MAX_angle_csx_hbond` parameter

**Code**

```python
pdb_file = 'hbond_S_2vh3_wConnect.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = sse_hydrogen_chalcogen_bond(traj, 0,1, MAX_angle_csx_hbond=150.0)
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'h-bond')
([], [[11, 5]])
([], [3.2067063450813293])
([], [[59.74052254307601, 3.1638450710847406, 149.59121339210552, 84.59585658972038]])
```

### 12.5. Error due to bad reconstruction of hydrogen positions

PDB2PQR is used to reconstruct hydogens in the PDb file. But in this case the H of the Cys is not reconstructed at the correct position. This leads to an error of analysis: the command return 'False' insted of 'h-bond'

**Code**

```python
pdb_file = 'hbond_S_5cpv_prq.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = sse_hydrogen_chalcogen_bond(traj, 0,1)
print(interaction.check_interaction)
print(interaction.get_atoms)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
False, False
([], [])
([], [])
([], [])
```



## 13. S/Se - Aromatic

**Code**

```python
pdb_file = 'aromatic_S_1azu.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = sse_aromatic(traj, 7,18) #HIS46 MET121
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'S-intermediate')
3.7583552070390196
51.534096128380675
```


## 14. Arg-Arg & Arg-aromatic

### 14.1. Arg-Aromatic - Parallel

**Code**

```python
pdb_file = 'ARG-ARG-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = arg_involved(traj, 0,1) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'Arg-Aromatic', 'parallel')
3.7319016658717015
6.318758521316121
```

### 14.2. Arg-Aromatic - Perpendicular

**Code**

```python
pdb_file = 'ARG-ARG-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = arg_involved(traj, 2,3) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'Arg-Aromatic', 'perpendicular')
4.710802378009887
83.06692941259753
```

### 14.3. Arg-Arg - Parallel

**Code**

```python
pdb_file = 'ARG-ARG-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = arg_involved(traj, 4,5) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'Arg-Arg', 'perpendicular')
4.070874452590942
72.61663724626918
```

### 14.4. Arg-Arg - Perpendicular

**Code**

```python
pdb_file = 'ARG-ARG-aromatic.pdb'
traj = md.load(pdb_file, top=pdb_file)
interaction = arg_involved(traj, 6,7) 
print(interaction.check_interaction)
print(interaction.get_distance)
print(interaction.get_angle)
```

**Result**

```
(True, 'Arg-Arg', 'parallel')
3.8860461115837097
17.299301309538308
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
