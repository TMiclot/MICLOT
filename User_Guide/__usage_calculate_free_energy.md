[User Guide home](Manual.md)

# Calculate free energy

## 1. Coulomb & Lennard-Jones

**coulomb_lj**(trajectory, res_index_A, res_index_B, force_field_path, *frame=0, coulomb_cutoff=12, lj_cutoff=12, lj_switch=8, solute_dielectric=1.0, solvent_dielectric=78.5*)

### Description

Calculate the Coulomb and Lennard-Jones (LJ) interaction between two residues.
It is possible to add a "hard" cutoff for both forces.
It is also possible calculate the Coulomb force using a reaction field, or to calculate the LJ with a switching function.

> [!CAUTION]
> It is important to have atom names consistent with the force field used: AMBER or CHARMM.
> The command integrates a simple conversion of atom names between these two force fields, but it may not be sufficient.
> Typically, an error due to atom names is:
> 
> AttributeError: 'NoneType' object has no attribute 'get'

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



## 2. Binding energy using the *contact based method*

**compute_binding_energy**(pdb_file_path, chainName_receptor, chainName_ligand, *temperature_celsius=25, write_outfile=True*)

### Description

Calculate the binding energy of a protein-protein complex using the [contact-based method]( https://doi.org/10.7554/eLife.07454).

### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| pdb_file_path       | string  | Path to the PDB file. | mandatory |
| chainName_receptor  | string  | List of chain names of the receptor. | mandatory |
| chainName_ligand    | string  | List of chain names of the ligand.   | mandatory |
| temperature_celsius | integer | Temperature in °C.          | optional |

### Returns

Return a Pandas Dataframe