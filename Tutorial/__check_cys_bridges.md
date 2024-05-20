[Tutorials home](Tutorials.md)

# Tutorial: How to check Cys-Cys bridges ?

This tutorial give you working example on how to check presence of cystein bridges into a PDB structure.

We will use the structure [5azz](https://www.rcsb.org/structure/5azz) of the seleno-insulin which contains both a diselenide bridge and two disulfide bridges.

## First step: load required modules

```python
import ***
```


## 1. Identify bridges between two Cys

> [!NOTE]
> Remember that energy calculation can be done only for disulfide bridge. The other cases use only distance parameters.

### 1.1. Disulfide bridge

**Code**

```python
pdb_file = "5azz.pdb"
traj = md.load(pdb_file, top=pdb_file)

bridge = cys_bridge(traj, 5, 10)

print(bridge.check_interaction)
print(bridge.get_energy)
print(bridge.get_energy_dse)
print(bridge.get_energy_dms)
print(bridge.get_distance)
print(bridge.get_angle)
```

**Result**

```
(True, 'disulfide')
(20.45071189975006, 3.8378589500717077)
((2.0735610679743317, 4.808189876333317), (7.254114775993074, 0.001970123592012989), 6.312876055857328)
((0.34683227647289283, 0.8042373286587984), 2.479890036667731, (0.1642325602320461, 0.042666748040239705))
(3.8848155736923218, 3.4517356753349304, 2.0521074533462524)
((-73.73775, -158.39505), (-105.7813, 59.413605), 62.37935, (111.06537, 116.4016))
```

### 1.2. Diselenide bridge

**Code**

```python
pdb_file = "5azz.pdb"
traj = md.load(pdb_file, top=pdb_file)

bridge = cys_bridge(traj, 6, 27)

print(bridge.check_interaction)
print(bridge.get_energy)
print(bridge.get_energy_dse)
print(bridge.get_energy_dms)
print(bridge.get_distance)
print(bridge.get_angle)
```

**Result**

```
(True, 'diselenide')
None
None
None
(4.719050824642181, 3.980471193790436, 1.878187507390976)
((171.5863, -60.29668), (64.167564, -87.63029), 90.05771, (115.554565, 112.92873))
```


### 1.3. NO bridge

**Code**

```python
pdb_file = "5azz.pdb"
traj = md.load(pdb_file, top=pdb_file)

bridge = cys_bridge(traj, 6, 19)

print(bridge.check_interaction)
print(bridge.get_energy)
print(bridge.get_energy_dse)
print(bridge.get_energy_dms)
print(bridge.get_distance)
print(bridge.get_angle)
```

**Result**

```
(False, False)
None
None
None
(17.158000469207764, 17.82523512840271, 18.588920831680298)
((171.5863, -57.02577), (115.526085, 13.11643), 76.61398, (115.554565, 116.53385))
```





## 2. Check the presence of Cys bridge on entire PDB file

### 2.1. Basic usage

**Code**

```python
cys_bridges_inPDB('5azz.pdb')
```

**Result**

```
Cys-Cys Bridge - Progress: 100%|█| 36/36 [00:00<0
```

Now you must have two new files *5azz_se-se_cys_bridges_inPDB.pdb* and *5azz_log.csv*. The first is the structure file in PDB format with renamed cysteines.
The other is the log file in CSV format. This table contains the sucessive pairs tested by the algorithm. When a pair is detected as bridging, the name of the two cysteines is changed accordingly: CYX or XSE.


### 2.2. Checking without writing logfile

**Code**

```python
cys_bridges_inPDB('5azz.pdb', logfile=False)
```

**Result**

```
Cys-Cys Bridge - Progress:  14%|▏| 5/36 [00:00<00

INFO - No logfile will be written.
	[+] 5 6 (False, False)
	[+] 5 39 (False, False)
	[+] 5 10 (True, 'disulfide')
	[+] 6 39 (False, False)
	[+] 6 19 (False, False)

Cys-Cys Bridge - Progress: 100%|█| 36/36 [00:00<0

	[+] 6 27 (True, 'diselenide')
	[+] 39 19 (True, 'disulfide')
```


### 2.3. Using without outputing a PDB file

**Code**

```python
cys_bridges_inPDB('5azz.pdb', outfile=False)
```

**Result**

```
Cys-Cys Bridge - Progress: 100%|█| 36/36 [00:00<0

INFO - No structure files (.pdb) have been modified or written.
```
