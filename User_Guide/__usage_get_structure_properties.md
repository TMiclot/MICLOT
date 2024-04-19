[User Guide home](Manual.md)

# Get structure properties

## 1. Sequence and Secondary structure

**get_sequence_secstruct**(pdb_file_path, *write_outfile=True*)

### Description

Return the sequence and secondary structure of each chain in a PDB.

The command returns two dictionaries:

- chainID and sequence
- chainID and secondary structure

The command also writes 2 output files (as CSV):

- residue names with their mdtraj toppology IDs and their secondary structure and their COM, COM of backbone and COM of side chain
- chainID and their corresponding sequence and secondary structure


> [!IMPORTANT]
> Depending on the PDB file, MDTraj may use different IDs for the same chain name in the PDB file.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| pdb_file_path | string  | Path to the PDB file.  | mandatory |
| write_outfile | boolean | If you want to write the output CSV files. | optional |

### Returns

The command return two dictionnaries: *dict_chainID_sequence, dict_chainID_ss*

| Return | Format | Description |
| ------ | ------ | --- |
| dict_chainID_sequence | dictionnary | chainID with coresponding sequence (string). |
| dict_chainID_ss | dictionnary | chainID with coresponding secondary structure (string). |





## 2. Protein region of amino acids *for a proteic complex*

