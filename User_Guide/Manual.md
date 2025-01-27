<img src="../__banner.png" alt="banner" class="center">

# User guide :notebook_with_decorative_cover:

## Installation

Please follow [installation instructions](__installation.md).



## Theory behind MICLOT

This section provides an in-depth explanation of the theoretical underpinnings behind the development of MICLOT. Understanding these concepts is essential for a comprehensive grasp of how to use the software.

- [Element properties](__element_properties.md)
- [Amino acids properties](__amino_acids_properties.md)
- [Cysteine bridges: Cys-Cys](__CysCys_bridges.md)
- [Non-bonding interactions](__nonbonding_interactions.md)
- [Free energy of residue pairs & Complex binding energy](__free_energy.md)
- [Glossary](__glossary.md)



## Usage

This section provides detailed information about the classes and methods available in MICLOT, along with guidance on how to utilize them effectively.

For a practical introduction to the software, you can also refer to the [tutorials](../Tutorial/Tutorials.md), which offer a hands-on exploration of MICLOT's capabilities.

- [Installation guide](__installation.md)
- [Prepare a structure](__usage_prepare_structure.md)
- [Get structure properties](__usage_get_structure_properties.md)
- [Identify non-bonded interaction](__usage_identify_nonbonded_interactions.md)
- [Identify Cys-Cys bridge](__usage_identify_CysCys_bridges.md)
- [Calculate free energy of interacting residues or of a complex](__usage_calculate_free_energy.md)
- [Data analysis](__usage_data_analysis.md) :warning: :construction:
- [Atoms & Amino acids Database](__usage_database.md)
- [Other tools](__usage_other_tools.md)



## Tips

### Reduce CPU usage

>[!TIP]
> For some cases you may want to limit the CPU usage, 2 ways are possible.
> 1. In Linux, you can restrain Python directy in Bash using `taskset`. For example, this command force python to be executed on CPU 0 to 3:
>
> ```bash
> taskset --cpu-list 0-3 python myscript.py
> ```
>
> 2. Directly in your Python script. For example you can allow Python using all CPU else the last two.
>
> ```python
> import os
> pid = os.getpid()
> max_cpu = os.cpu_count() -2 #Get total number of CPU -2 to avoid crash
> os.sched_setaffinity(pid, range(max_cpu)) 
> ```

### Kill all Python processes

>[!TIP]
> I case you over use your computer ressources, or want to kill the process, you can use this command:
>
> ```bash
> pkill -9 python
> ```



* * *
## Citing us
Miclot, T. & Timr, S. The famous title. *Journal* ... 
