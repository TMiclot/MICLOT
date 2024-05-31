<img src="../__banner.png" alt="banner" class="center">

# User guide :notebook_with_decorative_cover:

Description ... It's work so Well


## Introduction

This part provide a brief introduction to the software.

To start exploring the software, you can take look at the [tutorials](../Tutorial/Tutorials.md).


## Theory behind MICLOT

This part

- [Element properties](__element_properties.md)
- [Amino acids properties](__amino_acids_properties.md)
- [Cysteine bridges: Cys-Cys](__CysCys_bridges.md)
- [Canonical and non-canonical non-bonding interactions](__nonbonding_interactions.md)
- [Free energy of residue pairs & Complex binding energy](__free_energy.md)
- Surface shape analysis
- Network analysis
- [Glossary](__glossary.md)



## Usage

This part provide information about the classes and methods available in MICLOT.

- [Installation guide](__installation.md)
- [Prepare a structure](__usage_prepare_structure.md)
- [Get structure properties](__usage_get_structure_properties.md)
- [Identify non-bonded interaction](__usage_identify_nonbonded_interactions.md)
- [Identify Cys-Cys bridge](__usage_identify_CysCys_bridges.md)
- [Calculate free energy of interacting residues or of a complex](__usage_calculate_free_energy.md)
- [Atoms & Amino acids Database](__usage_database.md)
- [Other tools](__usage_other_tools.md)

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



* * *
## Citing us
Miclot, T. & Timr, S. The famous title. *Journal* ... 
