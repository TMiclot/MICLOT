[User Guide home](Manual.md)

# Identify Non-Bonded interactions

> [!NOTE]
> Each non-bonded interaction type is designed to be a [Classes](https://docs.python.org/3/tutorial/classes.html).
> For each classe it is possible to get somes properties and, for simplicity, keywords used are excalcy the same.



## x. Correlate PDB wth MDTraj topology

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





## x. Minimize structure using OpenMM

> [!NOTE]
> Minimization reporter come from [https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/report_minimization.html](https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/report_minimization.html)



### Arguments

> [!WARNING]
> The minimization algorithm don't constrain heavy atoms and a long iteration can lead to structure deformation. 