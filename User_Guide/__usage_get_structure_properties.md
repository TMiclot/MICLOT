[User Guide home](Manual.md)

# Get structure properties

All commands come from `miclot.utilities`

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

**get_protein_region**(pdb_file_path, chainID_receptor, chainID_ligand, *write_outfile=True*)

### Description

A command to return the protein region of all amino acids involved in a protein complex.
It use the maximum ASA reference values from: 

    - Thien et al. 2013   (https://doi.org/10.1371/journal.pone.0080635): Theoric
    - ibid. : Empiric
    - Miller et al. 1987  (https://doi.org/10.1016/0022-2836(87)90038-6)
    - Rose et al. 1985    (https://doi.org/10.1126/science.4023714)
    - Lins et al. 2003    (https://doi.org/10.1110/ps.0304803)
    - Samanta et al. 2002 (https://doi.org/10.1093/protein/15.8.659): Gly-X-Gly
    - ibid. : Ala-X-Ala
    - NACCESS software    (http://www.bioinf.manchester.ac.uk/naccess/)

This command return 2 dictionnaries:

    - ASA values for the complex, receptor and ligand. And the area of the interface. All values are in A^2.
    - chainID, of ligand and receptor, with their corresponding protein region for all amino acids (in "sequence" format). 

It also write 3 output files in CSV format:

    - Same as dictionnary: ASA values for the complex, receptor and ligand. And the area of the interface.
    - Same as dictionnary: chainID, of ligand and receptor, with their corresponding protein region for all amino acids (in "sequence" format).
    - Details information for each residues: SASA in complex and when free, rASA and drASA for each reference, then the protein region found.

The interface is calculated as follow:

$$
I = \frac {( {SASA}_{receptor} + {SASA}_{ligand} ) - {SASA}_{complex} }{2}
$$

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| pdb_file_path | string  | Path to the PDB file.  | mandatory |
| chainID_receptor | list of integer | MDTraj chain IDs of the receptor. | mandatory |
| chainID_ligand   | list of integer | MDTraj chain IDs of the ligand. | mandatory |
| write_outfile | boolean | If you want to write the output CSV files. | optional |

### Returns

The command return two dictionnaries: *dict_ASA_total_interface, dict_sequence_protein_region*

| Return | Format | Description |
| ------ | ------ | --- |
| dict_ASA_total_interface     | dictionnary | ASA values for the complex, receptor and ligand. And the area of the interface. |
| dict_sequence_protein_region | dictionnary | chainID, of ligand and receptor, with their corresponding protein region for all amino acids. |