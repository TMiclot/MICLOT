[Tutorials home](Tutorials.md)
# How to identify a given non-bonded interaction into an amino acid pair

## 1. C-bond
### Code
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

### Result
```
True
Atoms: [{1953, 1954, 1946, 1943}]
Angle: [[[169.06029, 170.11108], {1953, 1954, 1946, 1943}]]
Distance: [[2.871696650981903, {1953, 1954, 1946, 1943}]]
Energy: [[-21.651590464607594, {1953, 1954, 1946, 1943}]]
```
