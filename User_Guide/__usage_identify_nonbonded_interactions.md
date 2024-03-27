[User Guide home](Manual.md)

# Identify Non-Bonded interactions

> [!NOTE]
> Each non-bonded interaction type is designed to be a [Classes](https://docs.python.org/3/tutorial/classes.html).
> For each classe it is possible to get somes properties and, for simplicity, keywords used are excalcy the same.






# Identify Non-Bonded interactions between two amino acids

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
| .check_interaction | Check if the given interaction type exisit or not between the two amino acids. | Boolean (True / False) |  |
| .get_distance      | Return a list of distance of the C-bond, it also return the corresponding C-bond atom indices. | List of list: <br/> [[[distance],[indices]],...] | Å |
| .get_angle         | Return a list of angles of the C-bond, it also return the corresponding C-bond atom indices. <br/> $\theta 1$ angle: C...O=C <br/> $\theta 2$ angle: Z-C...O | Integer | degree |
| .get_energy        | Return a list of energy of the C-bond, it also return the corresponding C-bond atom indices. | List of list: <br/> [[[energy],[indices]],...]  | kJ/mol |
| .get_atoms         | Return a list of atoms index involved in C-bond. | List of set: <br/> [{indices},{indices},...] |  |





## 2. C5 Hydrogen bond - *inside the backbone of one residue*

**C5_hydrogen_bond**(trajectory, res_index, frame=0, *angular_tolerance=10.0, MAX_distance=2.7*)

### Description

Identify if their is a C5 H-bond interaction exist between the O and the H (or NH in CHARMM) lacated in the backbone of a residue.

### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| trajectory        | integer | MDTraj trajectory.  | mandatory |
| res_index         | integer | Index of the residue in MDTraj topology. | mandatory |
| angular_tolerance | integer | Set an absolute tolerance parameter N. So the range of angle is 140.0 +/- N. <br/> Default value of N: 10.0 <br/> See documentation concerning [numpy.isclose](https://numpy.org/doc/stable/reference/generated/numpy.isclose.html). | optional |
| frame             | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional  |
| MAX_distance      | integer | Maximum distance between O and H, in angstrom. <br/> Unit: Å <br/> Default value: 2.7 | optional  |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit. | Boolean (True / False) |  |
| .get_distance      | Distance between atoms O and H in the backbone of the residue. | distance | Å |
| .get_angle         | Values of Phi and Psi angle of the residue. | angle_phi, angle_psi | degree |





## 3. Hydrophobic interaction & Hydrophobe/Hydrophile clash

**hydrophobic**(trajectory, res_index_A, res_index_B, *frame=0, MAX_distance_COM=5.0, MAX_distance_CA=9.5, MIN_distance_CA=3.8*)

### Description

Interaction between tow hydrophobic residue or clash between hydrophilic/hydrophobic residues. Distance between their center of mass of the side chain must be $\leq$ 5.0 Å.

|             | Residues |
| ----------- | --- |
| Hydrophobic | ALA, CYS, ILE, LEU, MET, PHE, TRP, TYR, VAL |
| Hydrophilic | ARG, ASN, ASP, GLN, GLU, LYS |
> [!NOTE]
> Compare to [Pommié *et al.* (2004)](https://doi.org/10.1002/jmr.647) TYR is considered as hydrophobic, instead of neutral.


### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| trajectory  | integer | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional |
| MAX_distance_COM  | integer | Maximum distance between center of mass of each residue. <br/> Unit: Å <br/> Default value: 5.0 | optional  |
| MAX_distance_CA  | integer | Maximum distance between CA of each residue. <br/> Unit: Å <br/> Default value: 9.5 | optional  |
| MAX_distance_CA  | integer | Maximum distance between CA of each residue. <br/> Unit: Å <br/> Default value: 3.8 | optional  |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit. It return two boolean: the first determine if the interaction exist or not, the second identify *hydrophobic interaction* or *clash*. <br/> True,True: Exist & Hydrophobe <br/> True,False: Exist & Hydrophobic/Hydrophilic clash <br/> False,False: No interaction. | Boolean (True / False) |  |
| .get_distance      | Distances between CA-CA atoms and side chaine COM-COM of the two residues. | distance_CA, distance_COM | Å |
