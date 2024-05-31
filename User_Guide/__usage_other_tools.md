[User Guide home](Manual.md)

# Prepare a structure

All commands come from `miclot.utilities`


## 1. Correlate PDB with MDTraj topology

**mdtraj_chainID_2_chainName**(pdb_file_path, *write_outfile=True*)

### Description

Correlates a PDB file with its topology in MDTraj.

The command returns two dictionaries:

- convert chainID to chainName
- convert chainName to chainID

The command also writes three output files (as CSV):

- correspondence between mdtraj toppology and PDB file
- chainID and their corresponding chainName
- chainName and their corresponding chainID

> [!IMPORTANT]
> Depending on the PDB file, MDTraj may use different IDs for the same chain name in the PDB file.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| pdb_file_path | string  | Path to the PDB file.  | mandatory |
| write_outfile | boolean | If you want to write the output CSV files. | optional |

### Returns

The command return two dictionnaries: *dict_chainID_2_chainName, dict_chainName_2_chainID*

| Return | Format | Description |
| ------ | ------ | --- |
| dict_chainID_2_chainName | dictionnary | Convert MDTraj chainID to PDB chainName. |
| dict_chainName_2_chainID | dictionnary | Convert PDB chainName to MDTraj chainID. |