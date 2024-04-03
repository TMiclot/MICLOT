[User Guide home](Manual.md)
# Identify Cys-Cys bridges

To meet the user's needs, two commands have been implemented to identify cysteine weights:

1. The first option is to use a class that checks the formation of a bridge between two Cys and gives its type.
2. The second possibility is a command that checks over an entire PDB structure and writes a structure containing the names of the Cys modified if they make a bridge, with an output file (CSV format) that gives a general report on the bridge identification process.



## 1. Checking bridges between two Cys

**CysCys_bridge**(trajectory, trajectory, residue_index_A, residue_index_B, *frame=0, MAX_distance_XX=None, MAX_distance_CA=7.5, MIN_distance_CA=3.0, MAX_distance_CB=5.5, method=None, MAX_energy=75.0, MIN_energy=0.0*)

### Description



### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| trajectory  | integer | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |


### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ) |  |
| .get_distance      | Distances between CA-CA atoms of the two residues. | interger | Å |

> [!WARNING]
> This this method return same value for angle and energies as [Disulfide Bond Dihedral Angle Energy Server](https://services.mbi.ucla.edu/disulfide/).
> But it can return diffrent result compare to [Disulfide by Design 2.0 Server](http://cptweb.cpt.wayne.edu/DbD2/), because the angle measurements can be different between MDTraj, or a manual measurement with [VMD](https://www.ks.uiuc.edu/Research/vmd/), and what is given by this server.




## 2. Checking Cys-Cys bridges over a PDB structure

**check_cys_bridges**(pdb_file, *outfile=True, logfile=True*)

### Description



### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| pdb_file | string  | MDTraj trajectory. | mandatory |
| outfile | boolean | Index of residue A in MDTraj topology. | optional |
| logfile | boolean | Index of residue B in MDTraj topology. | optional |
