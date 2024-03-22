[User Guide home](Manual.md)
# Identify Non-Bonded interactions between two amino acids

Each non-bonded interaction type is designed to be a [Classes](https://docs.python.org/3/tutorial/classes.html).
For each classe it is possible to get somes properties and, for simplicity, keywords used are excalcy the same.


## 1. C-bond

**C_bond**(trajectory, res_index_A, res_index_B, *frame=0*)

### Description
Identify if their is C-bond interaction between C(sp3) with N or O or S.

### Arguments
| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| trajectory  | integer | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional  |

### Properties
| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit or not between the two amino acids. | Boolean: <br/> True / False | None |
| .get_distance      | Return a list of distance of the C-bond, it also return the corresponding C-bond atom indices. | List of list: <br/> [[[distance],[indices]],...] | nm |
| .get_angle         | Return a list of angles of the C-bond, it also return the corresponding C-bond atom indices. <br/> $\theta 1$ angle: C...O=C <br/> $\theta 2$ angle: Z-C...O | Integer | degree |
| .get_energy        | Return a list of energy of the C-bond, it also return the corresponding C-bond atom indices. | List of list: <br/> [[[energy],[indices]],...]  | kJ/mol |
| .get_atoms         | Return a list of atoms index involved in C-bond. | List of set: <br/> [{indices},{indices},...] | None |
