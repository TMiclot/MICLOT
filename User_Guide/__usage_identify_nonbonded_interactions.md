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
| .get_distance      | Return a list of distance of the C-bond, it also return the corresponding C-bond atom indices. | List of list: <br/> [[[distance],[indices]],...] < br/> distance and indices are integer. | Å |
| .get_angle         | Return a list of angles of the C-bond, it also return the corresponding C-bond atom indices. <br/> $\theta 1$ angle: C...O=C <br/> $\theta 2$ angle: Z-C...O | integer | degree |
| .get_energy        | Return a list of energy of the C-bond, it also return the corresponding C-bond atom indices. | List of list: <br/> [[[energy],[indices]],...] < br/> energy and indices are integer. | kJ/mol |
| .get_atoms         | Return a list of atoms index involved in C-bond. | List of set: <br/> [{indices},{indices},...] <br/> indices are integer. |  |





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
| .get_distance      | Distance between atoms O and H in the backbone of the residue. | integer | Å |
| .get_angle         | Values of Phi and Psi angle of the residue. | angle_phi, angle_psi | integer |





## 3. Hydrophobic interaction & Hydrophobe/Hydrophile clash

**hydrophobic**(trajectory, res_index_A, res_index_B, *frame=0, MAX_distance_COM=5.0, MAX_distance_CA=9.5, MIN_distance_CA=3.8*)

### Description

Identify interaction between tow hydrophobic residue or clash between hydrophilic/hydrophobic residues. Distance between their center of mass of the side chain must be $\leq$ 5.0 Å.

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
| MIN_distance_CA  | integer | Minimum distance between CA of each residue to consider as clash. <br/> Unit: Å <br/> Default value: 3.8 | optional  |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit. It return two boolean: the first determine if the interaction exist or not, the second identify *hydrophobic interaction* or *clash*. Output: <br/> True,True,"hydrophobic" <br/> True,False,"clash" <br/> False,False,None | Boolean (True / False / None) , String|  |
| .get_distance      | Distances between CA-CA atoms and side chaine COM-COM of the two residues. | distance_CA (integer), distance_COM (integer) | Å |





## 4. Charge clash & Charge repulsion *****

**charge_clash_repulsion**(trajectory, res_index_A, res_index_B, *frame=0, MAX_distance_CA=13.0, MIN_distance_charges=5.0*)

### Description

Check if their is a clash or a repulsion interaction between the two residues.

> [!NOTE]
> In cas of Arg-Arg stacking, this command will return a *clash* interaction. To have more presition concering Arg-Arg stacking you must use the dedicated command line.

| Residue | Charge | Charge location (*Not charged atoms*) |
| ------- | ------ | --- |
| LYS     | +      | NZ |
| ARG     | +      | CZ |
| HIP     | +      | ND1 |
| HSP     | +      | ND1 |
| ASP     | -      | CG |
| GLU     | -      | CD |

### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| trajectory  | integer | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional |
| MAX_distance_CA  | integer | Maximum distance between CA of each residue. <br/> Unit: Å <br/> Default value: 13.0 | optional  |
| MIN_distance_charges  | integer | Minimum distance between charges of each residue to consider as clash. <br/> Unit: Å <br/> Default value: 5.0 | optional  |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit. It return two boolean: the first determine if the interaction exist or not, the second identify *repulsion* or *clash*. Output: <br/> True,True,"clash" <br/> True,False,"repulsion" <br/> False,False,None | Boolean (True / False / None) , String|  |
| .get_distance      | Distances between CA-CA and charge-charge atoms of the two residues. | distance_CA (integer), distance_charge (,integer) | Å |






## 5. Salt bridge

**salt_bridge**(trajectory, res_index_A, res_index_B, *method="baker_hubbard", frame=0, MAX_distance_CA=13.0, MIN_distance_charges=4.0, distance_cutoff=3.0, angle_cutoff=120*)

### Description
Check if their is a strong electrostatic interaction involving an H-bond and an ionic bond between amino acids with opposite charges.
- Electrostatic interaction is given by distances.
- H-bonds can be computed using 3 diffrent methods: [baker_hubbard](https://mdtraj.org/1.9.4/api/generated/mdtraj.baker_hubbard.html), [kabsch_sander](https://mdtraj.org/1.9.4/api/generated/mdtraj.kabsch_sander.html), [wernet_nilsson](https://mdtraj.org/1.9.4/api/generated/mdtraj.wernet_nilsson.html).


| Residue | Charge | Charge location (*Not charged atoms*) | Atoms able to perform H-bond |
| ------- | ------ | --- | --- |
| LYS     | +      | NZ  | NZ HZ1 HZ2 HZ3 |
| ARG     | +      | CZ  | HE NE NH1 HH11 HH12 NH2 HH21 HH22 |
| HIP     | +      | ND1 | ND1 NE2 HE2 HD1 |
| HSP     | +      | ND1 | ND1 NE2 HE2 HD1 |
| ASP     | -      | CG  | OD1 OD2 HD2 |
| GLU     | -      | CD  | OE1 OE2 HE2 |

### Arguments

| Argument | Description | Format | Requirement |
| -------- | --- | --- | --- |
| trajectory  | integer | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Have no effect with the method 'baker_hubbard', because it implementation analyse the whole trajectory. <br/> Default value: 0 | optional |
| method      | string | Use 'baker_hubbard' or 'kabsch_sander' or 'wernet_nilsson' method to detect H-bond. | optional |
| MAX_distance_CA  | integer | Maximum distance between CA of each residue. <br/> Unit: Å <br/> Default value: 13.0 | optional  |
| MIN_distance_charges  | integer | Minimum distance between charges of each residue to consider as interacting. <br/> Unit: Å <br/> Default value: 4.0 | optional  |
| distance_cutoff | integer | In the 'baker_hubbard' method, the distance cutoff of Donor-H...Acceptor contact. <br/> Unit: Å <br/> Default value: 3.0 | optional |
| angle_cutoff    | integer | In the 'baker_hubbard' method, the angle cutoff of the angle $\theta$. <br/> Unit: degree <br/> Default value: 120 | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit. It return two boolean: the first determine if the ionic bond exist or not, the second identify hydrogen bond. Output: <br/> True,True: Interaction exist. <br/> True,False: Only ionic bond. <br/> False,True: Only H-bond. <br/> False,False: Interaction don't exist. | Boolean (True / False / None) & String|  |
| .get_distance      | Distances between CA-CA and charge-charge atoms of the two residues. | distance_CA (integer), distance_charge (integer) | Å |
| .get_hbond         | Return the atom indices forming H-bond in the pair. For each method, the output is converted to an **np.array**. | np.array |  |

> [!IMPORTANT]
> If you structure or trajectory don't contain hydrogens, the analysis can't provide H-bond analysis. So if their is a possible salt bride `.check_interaction` will return *True, False*. Of course this result can also depand the chosen method.


