[User Guide home](Manual.md)

# Calculate free energy

## Coulomb & Lennard-Jones

> [!IMPORTANT]  
> Here we use the files provided by [OpenMM](https://openmm.org/), because the data are presented in the same format for each force field.
> The command give the same values as [OpenMM](https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Computing%20Interaction%20Energies.html). Differences are lost in decimals.

> [!CAUTION]
> It is important to have atom names consistent with the force field used: AMBER or CHARMM.
> The command integrates a simple conversion of atom names between these two force fields, but it may not be sufficient.
> Typically, an error due to atom names is:
> 
> AttributeError: 'NoneType' object has no attribute 'get'


**coulomb_lj**(trajectory, res_index_A, res_index_B, force_field_path, *frame=0, coulomb_cutoff=12, lj_cutoff=12, lj_switch=8, solute_dielectric=1.0, solvent_dielectric=78.5*)

### Description

Calculate the Coulomb and Lennard-Jones interaction between two residues.
It is possible to add a "hard" cutoff for both forces.
It is also possible calculate the Coulomb force using a reaction field, or to calculate the L-J with a switching function.

### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| trajectory         | integer | MDTraj trajectory.  | mandatory |
| res_index_A        | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B        | integer | Index of residue B in MDTraj topology. | mandatory |
| force_field_path   | string  | Path to the force fiels AMBER or CHARMM. | mandatory |
| frame              | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| coulomb_cutoff     | integer | Cutoff applied to coulombic energy. It is used for *hard* cutoff and for reaction field. <br/> Unit: Å <br/> Default value: 12 | optional | 
| lj_cutoff          | integer | Cutoff applied to LJ energy. It is used for hard cutoff and for switching function. <br/> Unit: Å <br/> Default value: 12 | optional | 
| lj_switch          | integer |  <br/> Unit: Å <br/> Default value: 8 | optional |
| solute_dielectric  | integer |  <br/> Unit: Å <br/> Default value: 1.0 | optional |
| solvent_dielectric | integer |  <br/> Unit: Å <br/> Default value: 78.5 | optional |

> [!TIP]
> Because of their different file sizes, using AMBER ff14SB is faster than CHARMM36.
> Calculation with ff14SB run for 0.05 s while this time is 0.4 s for CHARMM.
> *These are approximate values which may vary depending on your computer configuration and the number of atoms in your selection.*

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ) |  |
| .get_distance      | Distances between CA-CA atoms of the two residues. | interger | Å |


