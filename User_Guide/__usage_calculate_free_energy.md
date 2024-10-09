[User Guide home](Manual.md)

# Calculate free energy

## 1. Coulomb & Lennard-Jones

All commands come from `miclot.coulomb_lj`

### 1.1. Without openMM

**coulomb_lj**(trajectory, res_index_A, res_index_B, force_field_path, *frame=0, coulomb_cutoff=12, lj_cutoff=12, lj_switch=8, solute_dielectric=1.0, solvent_dielectric=78.5*)

### Description

Calculate the Coulomb and Lennard-Jones (LJ) interaction between two residues.
It is possible to add a "hard" cutoff for both forces.
It is also possible calculate the Coulomb force using a reaction field, or to calculate the LJ with a switching function.

> [!CAUTION]
> It is important to have atom and residues names are consistent with the force field used: AMBER or CHARMM.
> The command integrates a simple conversion of atom names between these two force fields, but it may not be sufficient.
> Typically, an error due to atom names is:
> 
> AttributeError: 'NoneType' object has no attribute 'get'

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory         | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A        | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B        | integer | Index of residue B in MDTraj topology. | mandatory |
| force_field_path   | string  | Path to the force fiels AMBER or CHARMM. | mandatory |
| frame              | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| coulomb_cutoff     | integer | Cutoff applied to coulombic energy. It is used for *hard* cutoff and for reaction field. <br/> Unit: Å <br/> Default value: 12 | optional | 
| lj_cutoff          | integer | Cutoff applied to LJ energy. It is used for hard cutoff and for switching function. <br/> Unit: Å <br/> Default value: 12 | optional | 
| lj_switch          | integer | Switch distance applied to LJ energy. <br/> Unit: Å <br/> Default value: 8 | optional |
| solute_dielectric  | integer | Solute dielectric constant. <br/> Unit: Å <br/> Default value: 1.0 | optional |
| solvent_dielectric | integer | Solvent dielectric constant. <br/> Unit: Å <br/> Default value: 78.5 | optional |

> [!TIP]
> Because of their different file sizes, using AMBER ff14SB is faster than CHARMM36.
> Calculation with ff14SB run for 0.05 s while this time is 0.4 s for CHARMM.
> *These are approximate values which may vary depending on your computer configuration and the number of atoms in your selection.*

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .get_energy         | Return the total energy of the pair (LJ + Coulomb) without considering cutoff and with the "hard" cutoff for LJ and Coulomb. | integers: total, total_cutoff | kJ/mol |
| .get_energy_LJ      | Return all calculated energy for LJ. | integers: energy_LJ, energy_LJ_cutoff, energy_LJ_switch | kJ/mol |
| .get_energy_coulomb | Return all calculated energy for Coulomb. | integers: energy_C, energy_C_cutoff, energy_C_reactionfield | kJ/mol |

> [!TIP]
> The commands give the same values as [OpenMM](https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Computing%20Interaction%20Energies.html). Differences are lost in decimals.


### 1.2. Using openMM

**omm_coulomb_lj**(trajectory, index_residue_A, index_residue_B, *method=NoCutoff, nonbonded_cutoff=10.0, frame=0*)

### Description

>[!IMPORTANT]
>Adapted from https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/Computing%20Interaction%20Energies.html

Calculate Coulomb and Lennard-Jones energies using openMM. It automaticaly return values for AMBER and CHARMM force fields.

Sometime using this function avoid some error due to residues/atoms encountered by previous command `coulomb_lj`, but is much slower.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory      | mdtraj | MDTraj trajectory.  | mandatory |
| index_residue_A | integer | Index of residue A in MDTraj topology. | mandatory |
| index_residue_B | integer | Index of residue B in MDTraj topology. | mandatory |
| method          | openmm keyword  | Method used to calculate the energy. Allowed values: NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, or LJPME <br/> Default value: NoCutoff | optional |
| nonbonded_cutoff | float | Cutoff applied to nonbonded interactions.  <br/> Unit: Å <br/> Default value: 10 | optional |
| frame           | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .get_energy         | Return the total energy of the pair (LJ + Coulomb) without considering cutoff and with the "hard" cutoff for LJ and Coulomb. | float: amber_total, charmm_total | kJ/mol |
| .get_energy_LJ      | Return all calculated energy for LJ. | float: amber_lj, charmm_lj | kJ/mol |
| .get_energy_coulomb | Return all calculated energy for Coulomb. | float: amber_coulomb, charmm_coulomb | kJ/mol |




## 2. Binding energy using the *contact based method*

All commands come from `miclot.complex_binding`

**compute_binding_energy**(pdb_file_path, chainName_receptor, chainName_ligand, *temperature_celsius=25, write_outfile=True*)

### Description

Calculate the binding energy of a protein-protein complex using the [contact-based method](__free_energy.md#2-protein-binding-contacts-based-method). Currently the implementation allow user to work only with PDB ().
Because the method has been calibrated for use on fixed structures, but no mention has been made of molecular dynamics trajectories. So it was decided to restrict use to fixed structures: PDb, mmCIF, and all MDTraj readable files.

> [!IMPORTANT]
> It's important that the atom names match those of FreeSASA.
> If you analyse a structure extracted from a MD simulation, the atom names may not be well recognized. [PDBFixer](https://github.com/openmm/pdbfixer) is a good tool for this task.


### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| pdb_file_path       | string  | Path to the PDB file. | mandatory |
| chainName_receptor  | string  | List of chain names of the receptor. | mandatory |
| chainName_ligand    | string  | List of chain names of the ligand.   | mandatory |
| temperature_celsius | integer | Temperature in °C.          | optional |

### Returns

Return a Pandas Dataframe.