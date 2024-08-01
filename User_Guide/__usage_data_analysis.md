[User Guide home](Manual.md)

# Data analysis

All commands come from `miclot.analysis`


## 1. Correlate PDB with MDTraj topology

**plot_minimization**(path, *pdb_name=None, fig_size=[15, 4], linewidth=2, save_graph=False, dpi=300*)

### Description

Read the CSV file written durin the minimisation and plot:

- System energy
- Harmonic restraints energy
- Restraint force constant
- Constraint maximum relative error

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| pdb_file_path | string  | Path to the PDB file.  | mandatory |
| write_outfile | boolean | If you want to write the output CSV files. | optional |

### Returns

The command return the graph in Matplotlib format.