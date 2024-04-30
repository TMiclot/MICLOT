[User Guide home](Manual.md)

# Prepare a structure

A set of functions present in the *utilities* used to prepare a structure for analysis.


## x. Add missing atoms & renam residues by their protonated state
**pbd2pqr_parse**(pdb_file_path, *force_field='AMBER', ph=7.0, write_logfile=True*)

### Description

It is use to prepare structure for analysis. It is able to add missing atoms and find protonation state of redidues.

The function is a parser to PDB2PQR. It return the corresponding PQR file of the structure and a PDB file corresponding to the PQR file
(thank to the option --pdb-output). Future analysis must be done on this PDB output file.

Option used are: titration method is [PROPKA](https://github.com/jensengroup/propka), it keep chain IDs in the PQR file and remove water molecules.

> [!NOTE]
> The command use 'subsystem' insted of the Python API, because as mentioned in the documentation:
> > The API is still changing and there is currently no guarantee that it will remain stable between minor releases.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| pdb_file_path | string  | Path to the PDB file. | mandatory |
| write_logfile | boolean | If you want to write the log of PDB2PQR into a text file. | optional |
| force_field   | string  | Force field to use 'AMBER' or 'CHARMM'. <br/> Default value: 'AMBER' | optional |
| ph            | float   | pH value use to identify protonation state. <br/> Default value: 7.0 | optional |

### Returns

The command return 3 files:

| Name   | Description |
| ------ | --- |
| *PDBname*_pdb2pqr.log | Log file of PDB2PQR |
| *PDBname*_pqr.pdb | Output PDB file with the protonated state names. |
| *PDBname*.pqr     | Output PRQ file. |



## x. Correlate PDB with MDTraj topology

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

### Description

> [!NOTE]
> Minimization reporter come from [https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/report_minimization.html](https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/report_minimization.html)

This command is a parser to minimize a structure using [openMM](https://openmm.org/). It generate a PDB file and an output file in CSV format containing the energy of the structure at each step of the minimization:

- system energy: the current potential energy of the system
- restraint energy: the energy of the harmonic restraints
- restraint strength: the force constant of the restraints (in kJ/mol/nm^2)
- max constraint error: the maximum relative error in the length of any constraint

For more details, the options: `nonbondedMethod=NoCutoff, constraints=HBonds` are used. The first ensure all Coulomb and Lennar-Jones force a re calculated without cutoff, the second ensure *the lengths of all bonds that involve a hydrogen atom are constrained*.

> [!WARNING]
> The minimization algorithm don't constrain heavy atoms and a long minimization can lead to structure deformation.

### Arguments 

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| pdb_file_path | string  | Path to the PDB file. | mandatory |
| ff | string | Force field to use. Values are 'amber' or 'charmm'. <br/> Default value: 'amber' | optional |
| max_iterations | integer | Maximum number of iterations. <br/> Default value: 100 | optional |

### Returns

The command return a minimization log file, in CSV format, and the minimized structure in PDB format.