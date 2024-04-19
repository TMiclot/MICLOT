[Tutorials home](Tutorials.md)

# Get structure properties from a PDB file

## 1. Sequence and Secondary structure

**Code**

```python
chainID_sequence, chainID_ss = get_sequence_secstruct('5azz.pdb')
print("chainID_sequence:", chainID_sequence)
print("chainID_ss:", chainID_ss)
```

**Result**

```
chainID_sequence: {0: 'GIVEQCUASVCSLYQLENYCN', 1: 'FVNQHLUGSHLVEALYLVCGERGFFYTPK', 2: 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX', 3: 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'}
chainID_ss: {0: 'CHHHHHCCCCCCHHHHHCCEC', 1: 'CCCCCCCCHHHHHHHHHHHHHHCECCCCC', 2: 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN', 3: 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'}
```

The command also output 2 CSV files:

- 5azz_residues_secondaryStructure_and_COMs.csv
- 5azz_sequence_secondaryStructure.csv





## 2. Protein region of amino acids *for a proteic complex*