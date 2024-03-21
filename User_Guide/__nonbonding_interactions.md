[User Guide home](Manual.md)
# Identify Non-Bonded interactions between two amino acids

Each non-bonded interaction type is designed to be a [Classes](https://docs.python.org/3/tutorial/classes.html).
For each classe it is possible to get somes properties and, for simplicity, keywords used are excalcy the same.
Else for some interaction classes which requires special properties.


## 1. C-bond

**C_bond**(trajectory, res_index_A, res_index_B, *frame=0*)

### Description
Interaction between C(sp3) with N or O or S.

### Arguments
| Argument | Description | Format | Required |
| -------- | --- | --- | --- |
| trajectory  | integer | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional  |

### Properties
| Property | Effect | Return |
| -------- | --- | --- |
| .check_interaction | Check if the given interaction type exisit or not between the two amino acids. |  |
| .get_distance      |  |  |
| .get_angle         |  |  |
| .get_energy        |  |  |
| .get_atoms         |  |  |
