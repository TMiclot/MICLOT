[User Guide home](Manual.md)
# Identify Cys-Cys bridges

To meet the user's needs, two commands have been implemented to identify cysteine weights:

1. The first option is to use a class that checks the formation of a bridge between two Cys and gives its type.
2. The second possibility is a command that checks over an entire PDB structure and writes a structure containing the names of the Cys modified if they make a bridge, with an output file (CSV format) that gives a general report on the bridge identification process.



## 1. Checking bridges between two Cys

**CysCys_bridge**(trajectory, trajectory, residue_index_A, residue_index_B, *frame=0, MAX_distance_XX=None, MAX_distance_CA=7.5, MIN_distance_CA=3.0, MAX_distance_CB=5.5, method=None, MAX_energy=75.0, MIN_energy=0.0*)

### Description

 Identify disulfide, diselenide or selenosulfide bridge between two Cys.
 It ca use two methods to calculate the bridge energy:

 - Dihedral strain energy (DSE)
 - Disulfide model structure (DMS)



### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| trajectory  | integer | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| MAX_distance_XX | integers. | Maximum distance beween S, Se, atoms. <br/> Unit: Å <br/> Default value: Covalent radii | optional |
| MAX_distance_CA | integer | Maximum distance between the CA of each Cys. <br/> Unit: Å <br/> Default value: 7.5 | optional |
| MIN_distance_CA | integer | Minimum distance between the CA of each Cys. <br/> Unit: Å <br/> Default value: 3.0 | optional |
| MAX_distance_CB | integer | Maximum distance between the CB of each Cys. <br/> Unit: Å <br/> Default value: 5.5 | optional |
| method | boolean | <br/> Default value: 0 | optional |
| MAX_energy | integer | Maximum energy to consider a bridge. To use it, remember to choose a method. <br/> Unit: kcal/mol or kJ/mol (Depend on method) <br/> Default value: 0 | optional |
| MIN_energy | integer | Minimum energy to consider a bridge. To use it, remember to choose a method. <br/> Unit: kcal/mol or kJ/mol (Depend on method)  <br/> Default value: 0 | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ) |  |
| .get_distance      | Return all CA-CA, CB-CB and X-X distances. | interger | Å |
| .get_angle         | Return all X1, X2, X3 and theta angles.    | interger | degree |
| .get_energy | Return the calculated energy for the method DSE and for the method DMS. | dse_total (integer), dms_total (integer) | DSE in kJ/mol, DMS  in kcal/mol |
| .get_energy_dse | Return the calculated energies for: X1 angles, X2 angles and X3 angle used in the method DSE. | integers: (X1_A,  X1_B), (X2_A,  X2_B), X3 | kJ/mol |
| .get_energy_dms | Return the calculated energies for: X1 angles, X2 angles and X3 angle used in the method DSE. | integers: (X1_A,  X1_B), X3, (theta_A,  theta_B) | kcal/mol |

> [!WARNING]
> This this method return same value for angle and energies as [Disulfide Bond Dihedral Angle Energy Server](https://services.mbi.ucla.edu/disulfide/).
> But it can return diffrent result compare to [Disulfide by Design 2.0 Server](http://cptweb.cpt.wayne.edu/DbD2/), because the angle measurements can be different between MDTraj, or a manual measurement with [VMD](https://www.ks.uiuc.edu/Research/vmd/), and what is given by this server.




## 2. Checking Cys-Cys bridges over a PDB structure

**check_cys_bridges**(pdb_file, *outfile=True, logfile=True*)

### Description

This command check a PDB structure to identify all CYS-CYS bridges. It can identify disulfide, diselenide or selenosulfide bridges.
(seleno-)cystein are renamed: CYS to CYX, or SEC to XSE if they are involved into a bridge.

> [!TIP]
> A progress bar is displayed to inform the user on the status of the analysis.

### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| pdb_file | string  | PDB structure path. | mandatory |
| outfile  | boolean | Write a PDB file with modified Cys names if any bridge is detected. <br/> Default value: True | optional |
| logfile  | boolean | Generate an output file in CSV format containing all tests performed and their results. If set to False the result will be print in the terminal, but with less informations. <br/> Default value: True | optional |

> [!IMPORTANT]
> In the output PDB file, chain name and residue number can be renamed and/or renumbered, because that's how MDTraj works.
> Take these changes into account in your analysis.

### Information contained in the `logfile`

The logfile is in [CSV format](https://en.wikipedia.org/wiki/Comma-separated_values), it contain 12 columns.
Part of them are related to the first cystein (_A columns), or to the second cystein (_B columns).

The logfile show all Cys-Cys pairs tested and their results. 
In this table, (seleno-)cystein are renamed only when they are involved into a bridge.

| Column name | Description |
| ----------- | ----------- |
| ResName     | Name of the cysteine. |
| ResID       | ID of the cysteine in the MDTraj topology. |
| ResSeq      | Residue number in the sequence. |
| ChainID     | ID of the residue chain in the MDTraj topology. |
| ChainName   | Converted ID number to alphabetic name. |
| Bridge      | A boolean value which is True if their is a bridge between the two cysteine, or which is False if not. |
| Type        | Type of the bridge: disulfide, diselenide, selenosulfide or Fals if their is no bridge. |



