[User Guide home](Manual.md)

# Identify Non-Bonded interactions

All commands come from `miclot.interactions`

> [!NOTE]
> Each non-bonded interaction type is designed to be a [Classes](https://docs.python.org/3/tutorial/classes.html).
> For each classe it is possible to get somes properties and, for simplicity, keywords used are excalcy the same.

> [!IMPORTANT]
> Keep in mind that protonation state name is not always used in a structure file. Depending on the context, it may be useful to use PDB2PQR to get this information.

> [!WARNING]
> 1. If you analyse PDB structure file, be sure to have the [CONNECT](https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html) section in the file, or the [struct_conn](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/struct_conn.html) category for PDBx/mmCIF files. 
> Alternatively you can use the [`fix_topology_bonds`](__usage_prepare_structure.md#2-create-a-trajectory-with-correct-bonds-in-the-topology) to generate a corrected topology.
> <br/> It is mandatory if your file contain not [Standard code](__amino_acids_properties.md#2-the-3-letter-codes-of-standard-residues-in-pdb) residue names. For example, bonds of seleno-methionine (MSE) are not recognized if the file don't containt the information about connectivity between atoms.
>
> 2. It is also necessary to have hydrogens atoms to identify hydrogen bonds, else for the S/Se mediated H-bond.



## Identify all interfactions in a protein, or a complex

**interaction_table_whole_system**(trajectory, *list_pairs="all", frame=0, MAX_distance_contact=14.0, use_tqdm=False, write_outfile=True, path_name_outfile="interaction_table_whole_system.csv"*):
    
### Description

*Here a protein or a protein complex is named system, without difference.*

In a system analyse all possible pairs (or given pairs) and return a complete table of their interaction types, with: 

- Coulomb and Lennard-Jones energyes as calculated by openMM using both AMBER and CHARMM,
- Presence/Absence of C5-Hbond for each residue perfoming the pair.

The command return a Pandas DataFrame with all data and write it as CSV file.


>[!TIP]
> You can easily create a list with all possible pairs between given residues, without redundancy, using [itertools](https://docs.python.org/3/library/itertools.html). Below, a very short example:
>
> ```python
> import itertools
> 
> list_residues = [1,2,3,4]
>
> # Create a list with all possible pairs between residues IDs: 1,2,3,4
> # (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)
> list_pairs = list(itertools.combinations(list_residues, 2))
> ```

> [!CAUTION]
> 1. This command use all CPUs seen by Python, so it can use all CPUs of your computer. To ensure avoid any crash, be sure to limit the number of CPUs used by: Jupyter, or your script. To do this, you can take look at this [Tip](Manual.md#usage)
> 2. This command can increase the memory usage by ~ 2/3 Go. Be sure your computer has enough free memory.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional  |
| list_pairs  | list | List of pair of residue. By default take all possible pair in the system. Else can be set by user, using this format: [[id1,id2], [id3,id4], ...] | optional |
| MAX_CA_distance | float | Maximum distance between CA atom of the two residue. Pairs with a greater distance will be ignored. <br/> Default value: 14.0 angstrom | optional |
| write_outfile | boolean | Write the finale interaction table in CSV format (True/False). <br/> Default value: True | optional |
| path_name_outfile | string | Path, containing the name, of the exported CSV table. By default the file is exported in the curent location and the name: "interaction_table_whole_system.csv" <br/> Custom: "my/path/interaction_table_whole_system.csv" | optional |
| use_tqdm | boolean | Display tqdm proggress bar (True/False). <br/> Default value: False | optional |


## Identify all interfactions for a pair of residues

**identify_all_interaction_pair**(trajectory, pair, *frame=0*)

### Description

A fucntion to identify all interaction between residues forming a pair.

It run all interaction commands (with default parameters) and return a Pandas DataFrame with informations of each residue, the number of interaction doing by the pair, a detail of the interaction type and subtype.
If the interaction exist: 1, else 0. If the interaction can't exist: NaN.

>[!IMPORTANT]
> Exception for **salt bridges**: it return 1 if both H-bond and short distance between charges are identified, or return 0.5 if only one of them is identified.


### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| pair        | list | Pair of residues. The format is a list of MDTraj indices. <br/> Example: [0,3] | mandatory |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional  |



## 1. C-bond

**C_bond**(trajectory, res_index_A, res_index_B, *frame=0, MAX_distance=3.6, MIN_distance=2.5, \
                 MAX_angle_COC=180.0, MIN_angle_COC=160.0, MAX_angle_ZOC=180.0, MIN_angle_ZOC=160.0, \
                 MAX_energy=-2.0, MIN_energy=-22.2*)

### Description

Identify if their is C-bond interaction between C(sp3) with O.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional  |
| MAX_distance  | float | Maximum C-O distance. <br/> Default value: 3.6 angstrom |
| MIN_distance  | float | Mimimum C-O distance. <br/> Default value: 2.5 angstrom |
| MAX_angle_COC | float | Maximum C..O=C angle. <br/> Default value: 180.0 degree |
| MIN_angle_COC | float | Minimum C..O=C angle. <br/> Default value: 160.0 degree |
| MAX_angle_ZOC | float | Maximum Z-O...C angle. <br/> Default value: 180.0 degree |
| MIN_angle_ZOC | float | Minimum Z-O...C angle. <br/> Default value: 160.0 degree |
| MAX_energy    | float | Maximum energy of the interaction. <br/> Default value: -2.0 kJ/mol |
| MIN_energy    | float | Minimum energy of the interaction. <br/> Default value: -22.0 kJ/mol |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit or not between the two amino acids. | Boolean (True / False) |  |
| .get_distance      | Return a list of distance of the C-bond, it also return the corresponding C-bond atom indices. | List of list: <br/> [[[distance],[indices]],...] < br/> distance and indices are float. | Å |
| .get_angle         | Return a list of angles of the C-bond, it also return the corresponding C-bond atom indices. <br/> $\theta 1$ angle: C...O=C <br/> $\theta 2$ angle: Z-C...O | float | degree |
| .get_energy        | Return a list of energy of the C-bond, it also return the corresponding C-bond atom indices. | List of list: <br/> [[[energy],[indices]],...] < br/> energy are float and indices are integer. | kJ/mol |
| .get_atoms         | Return a list of atoms index involved in C-bond. | List of set: <br/> [{indices},{indices},...] <br/> indices are integer. |  |





## 2. C5 Hydrogen bond - *inside the backbone of one residue*

**C5_hydrogen_bond**(trajectory, res_index, frame=0, *angular_tolerance=10.0, MAX_distance=2.7*)

### Description

Identify if their is a C5 H-bond interaction exist between the O and the H (or NH in CHARMM) lacated in the backbone of a residue.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory        | mdtraj | MDTraj trajectory.  | mandatory |
| res_index         | integer | Index of the residue in MDTraj topology. | mandatory |
| angular_tolerance | float | Set an absolute tolerance parameter N. So the range of angle is 140.0˚ +/- N. <br/> Default value of N: 10.0 <br/> See documentation concerning [numpy.isclose](https://numpy.org/doc/stable/reference/generated/numpy.isclose.html). | optional |
| frame             | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional  |
| MAX_distance      | float | Maximum distance between O and H, in angstrom. <br/> Unit: Å <br/> Default value: 2.7 | optional  |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit. | Boolean (True / False) |  |
| .get_distance      | Distance between atoms O and H in the backbone of the residue. | float | Å |
| .get_angle         | Values of Phi and Psi angle of the residue. | angle_phi, angle_psi | float |





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

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional |
| MAX_distance_COM  | float | Maximum distance between center of mass of each residue. <br/> Unit: Å <br/> Default value: 5.0 | optional  |
| MAX_distance_CA  | float | Maximum distance between CA of each residue. <br/> Unit: Å <br/> Default value: 9.5 | optional  |
| MIN_distance_CA  | float | Minimum distance between CA of each residue to consider as clash. <br/> Unit: Å <br/> Default value: 3.8 | optional  |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit. It return two boolean: the first determine if the interaction exist or not, the second identify *hydrophobic interaction* or *clash*. Output: <br/> True,True,"hydrophobic" <br/> True,False,"clash" <br/> False,False,None | Boolean (True / False / None) , String|  |
| .get_distance      | Distances between CA-CA atoms and side chaine COM-COM of the two residues. | distance_CA (float), distance_COM (float) | Å |





## 4. Charge clash & Charge repulsion

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

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Default value: 0 | optional |
| MAX_distance_CA  | float | Maximum distance between CA of each residue. <br/> Unit: Å <br/> Default value: 13.0 | optional  |
| MIN_distance_charges  | float | Minimum distance between charges of each residue to consider as clash. <br/> Unit: Å <br/> Default value: 5.0 | optional  |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit. It return two boolean: the first determine if the interaction exist or not, the second identify *repulsion* or *clash*. Output: <br/> True,True,"clash" <br/> True,False,"repulsion" <br/> False,False,None | Boolean (True / False / None) , String|  |
| .get_distance      | Distances between CA-CA and charge-charge atoms of the two residues. | distance_CA (float), distance_charge (,float) | Å |






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

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Have no effect with the method 'baker_hubbard', because it implementation analyse the whole trajectory. <br/> Default value: 0 | optional |
| method      | string | Use 'baker_hubbard' or 'kabsch_sander' or 'wernet_nilsson' method to detect H-bond. | optional |
| MAX_distance_CA  | float | Maximum distance between CA of each residue. <br/> Unit: Å <br/> Default value: 13.0 | optional  |
| MIN_distance_charges  | float | Minimum distance between charges of each residue to consider as interacting. <br/> Unit: Å <br/> Default value: 4.0 | optional  |
| distance_cutoff | float | In the **baker_hubbard** method, the distance cutoff of Donor-H...Acceptor contact. <br/> Unit: Å <br/> Default value: 3.0 | optional |
| angle_cutoff    | float | In the **baker_hubbard** method, the angle cutoff of the angle $\theta$. <br/> Unit: degree <br/> Default value: 120 | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit. It return two boolean: the first determine if the ionic bond exist or not, the second identify hydrogen bond. Output: <br/> True,True: Interaction exist. <br/> True,False: Only ionic bond. <br/> False,True: Only H-bond. <br/> False,False: Interaction don't exist. | Boolean (True / False / None) & String|  |
| .get_distance      | Distances between CA-CA and charge-charge atoms of the two residues. | distance_CA (float), distance_charge (float) | Å |
| .get_hbond         | Return the atom indices forming H-bond in the pair. For each method, the output is converted to an np.array. | np.array |  |

> [!IMPORTANT]
> If your structure, or trajectory, don't contain hydrogens, the analysis can't provide H-bond analysis. So if their is a possible salt bride `.check_interaction` will return *True, False*. Of course this result can also depand of the chosen method.



## 6. Hydrogen bond

**hydrogen_bond**(trajectory, res_index_A, res_index_B, *method="baker_hubbard", distance_cutoff=3.0, angle_cutoff=150, frame=0*)

### Description

Identify H-bond in a structure, or a trajectory. H-bonds can be computed using 3 diffrent methods: [baker_hubbard](https://www.mdtraj.org/1.9.8.dev0/api/generated/mdtraj.baker_hubbard.html), [kabsch_sander](https://www.mdtraj.org/1.9.8.dev0/api/generated/mdtraj.kabsch_sander.html), [wernet_nilsson](https://www.mdtraj.org/1.9.8.dev0/api/generated/mdtraj.wernet_nilsson.html).

Based on distances, this class can also discriminate regular, low-barrier, and single-well H-bond subtypes. This information is directly related to the energy profile of the interaction.

> [!IMPORTANT]
> - The structure, or a trajectory must contain hydrogens.
> - To avoid redundancy with the class [salt_bridge](#5-salt-bridge), this command exclude H-bond involved in a salt bridge only.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. </br> Have no effect with the method 'baker_hubbard', because it implementation analyse the whole trajectory. <br/> Default value: 0 | optional |
| method      | string | Use 'baker_hubbard' or 'kabsch_sander' or 'wernet_nilsson' method to detect H-bond. | optional |
| distance_cutoff | float | In the **baker_hubbard** method, the distance cutoff of Donor-H...Acceptor contact. <br/> Unit: Å <br/> Default value: 3.0 | optional |
| angle_cutoff    | float | In the **baker_hubbard** method, the angle cutoff of the angle $\theta$. <br/> Unit: degree <br/> Default value: 120 | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | ----------- | ------ | ---- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ) |  |
| .get_atoms         | Return a list of atoms index involved in H-bond. | numpy.array: <br/> [{indices},{indices},...] |  |
| .get_distance      | Return a list of distance between donnor and acceptor, it also return the corresponding C-bond atom indices.  | List: <br/> [[distance,[indices]],...] | Å |
| .get_angle         | Return a list of angles (DHA) of all H-bonds, it also return the corresponding atom indices. | List of float: <br/> [[angle,[indices]],...] | degree |
| .get_subtype       | Return the subtype of the H-bond.  | List: <br/> [[subtype,[indices]],...] | |



## 7. van der Waals

**van_der_waals**(trajectory, res_index_A, res_index_B, frame=0, set_hydrogen=True, distance_tolerance=0.5, MIN_contact_numbers=1)

### Description

Identify van der Waals interaction between two residues.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| set_hydrogen | boolean (True / False) | Set if hydrogens are taken in acount, or not, for vdW interaction. <br/> Default value: True | optional |
| distance_tolerance | float | The range of distance is $radii_{vdw \space atom \space 1} + radii_{vdw \space atom \space 2} + N$ <br/> Range of values: 0.0 to 0.6 Å. <br/> Default value: 0.5 Å | optional |
| MIN_contact_numbers | float | Min number of contacts tow residues must have to consider a vdW interaction. <br/> Default value: 1 | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction   | Check if the given interaction type exisit.  | Boolean (True / False ) & string() |  |
| .get_distance        | Return the list distances between atoms pairs making vdW interaction and the list of their index. | list_distance (list), list_contacts (list) | Å |
| .get_number_contacts | Return the number of atom-atom vdW contact between the two residues. | integer |  |
| .get_interface       | Return the interface contact between the two residues. *None* is return when the vdw interaction don't exist. The value is given by the equation: <br/> $SASA_{residu \space A} + SASA_{residu \space B} - SASA_{pair \space AB}$ | float or boolean | $Å^2$ |





## 8. Amino- $\pi$

**amino_pi**(trajectory, res_index_A, res_index_B, *frame=0, MAX_distance=5.5, angular_tolerance=30.0*)

### Description

Identify the interaction between the amino group of Asn or Gln and the $\pi$ ring of the residue.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory   | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A  | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B  | integer | Index of residue B in MDTraj topology. | mandatory |
| frame        | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| MAX_distance | float | Maximum distance between the atom N and the center of mass (COM) of the ring. <br/> Unit: Å <br/> Default value: 5.5 | optional |
| angular_tolerance | float | Set an absolute tolerance parameter N. So the range of angle is 90.0˚ +/- N. <br/> Default value of N: 30.0 <br/> See documentation concerning [numpy.isclose](https://numpy.org/doc/stable/reference/generated/numpy.isclose.html).  | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ) |  |
| .get_distance      | Distances between the COM of the ring and the N atom. | float | Å |
| .get_angle         | Angle between vector normal of the aromatic ring plan and the vector COM $\rightarrow$ N. | float | degree |






## 9. Charge-Aromatic ring

**charge_aromatic**(trajectory, res_index_A, res_index_B, frame=0, MAX_distance=5.5, MIN_pi_angle=60.0, MAX_quadrupole_angle=35.0)

### Description

Identify the interaction with the charged amino acid and an aromatic rings.
It identify 3 subtypes, where *charge* is cation or anion:

- *charge*-$\pi$
- *charge*-intermediate (when the angle is between the $\pi$ area and the quadrupole area)
- *charge*-quadrupole

>[!NOTE]
> Protonated histidine (HIP or HSP) are considered as charged residue only.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory   | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A  | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B  | integer | Index of residue B in MDTraj topology. | mandatory |
| frame        | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| MAX_distance | float | Maximum distance between the charge and the center of mass (COM) of the ring. <br/> Unit: Å <br/> Default value: 5.5 | optional |
| MIN_pi_angle | float | Minumum angle to set Pi area. (The max angle is 90) <br/> Unit: degree <br/> Default value: 60.0 | optional |
| MAX_quadrupole_angle | float | Maximum angle to set quadrupole area. (The min angle is 0) <br/> Unit: degree <br/> Default value: 35.0 | optional |


### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ) |  |
| .get_distance      | Distances between the COM of the ring and the charged atom. | float | Å |
| .get_angle         | Angle between vector normal of the aromatic ring plan and the vector COM $\rightarrow$ charge. | float | degree |




## 10. Aromatic-Aromatic

**aromatic_aromatic**(trajectory, res_index_A, res_index_B, *frame=0, MAX_angle_planarity=30.0, MAX_distance_COM=5.5, MIN_distance_offset=1.6, MAX_distance_offset=2.0, MIN_pi_angle=60.0, MAX_quadrupole_angle=35.0, MAX_angle_Tshaped=5.0*)

### Description

Identify the stacking of two aromatic residues. Their is parallel, offset, t-shaped, y-shaped, coplanar, and when the angular parameters don't correspond to a specific subtype: intermediate.

Please note that protonated histidine (HIP or HSP) are not taken in acount are not taken into account because they can make charge-aromatic interactions.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| MAX_distance_COM | float | Distance between the COMs of each residue <br/> Default value: 5.5 Å | optional |
| MIN_distance_offset | float | Minimum distance between the COM of the first residue and the COM of the second residue projected on the plane of the first residue. <br/> Default value: 1.6 Å | optional |
| MAX_distance_offset | float | Maximum distance between the COM of the first residue and the COM of the second residue projected on the plane of the first residue. <br/> Default value: 2.0 Å | optional |
| MIN_pi_angle | float | Minimum angle defining the Pi area (the maximum is 90˚). <br/> Default value: 60.0˚ | optional |
| MAX_quadrupole_angle | float | Maximum angle defining the quadrupole area (the maximum is 0˚). <br/> Default value: 35.0˚ | optional |
| MAX_angle_planarity | float | Maximum angle defining the planarity between the two planes (the maximum is 0˚). <br/> Default value: 30.0˚ | optional |
| MAX_angle_Tshaped | float | Maximum angle between the normal vector of an aromatic plane of residue 1 and the vector $\overrightarrow{COM \space C}$ of the ring of the residue 2 (the maximum is 0˚). <br/> Default value: 5.0˚ | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ), string (subtype) |  |
| .get_angle         | Angle between the two plane, angle beteewn the vector COM $\rightarrow$ COM and plane A then with plane B. Finally it return the hsaped position angle.| float | degree |
| .get_distance      | Return the distances between COM-COM and distance between the COM and the projected COM on a plane. | float | Å |




## 11. ARG involved

**arg_involved**(trajectory, res_index_A, res_index_B, *frame=0, MAX_distance=6.0, MIN_pi_angle=60.0, MAX_quadrupole_angle=35.0*)

### Description

A specific class to precise the ARG-Aromatic and ARG-ARG interaction. It is possible to identify 3 subtypes: perpendicular, parallel, intermediate.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| MAX_distance | float | Distance between the CZ atom of the ARG and COM of and aromatic ring or CZ of another ARG. <br/> Default value: 6.0 Å | optional |
| MIN_pi_angle | float | Minimum angle defining the Pi area (the maximum is 90˚). <br/> Default value: 60.0˚ | optional |
| MAX_quadrupole_angle | float | Maximum angle defining the quadrupole area (the maximum is 0˚). <br/> Default value: 35.0˚ | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ) |  |
| .get_angle         | Angle between the 2 Arg plane or the angle between the Arg plane and the aromatic plane. | float | degree |
| .get_distance      | Distance between the CZ of each Arg, of distance between the CZ of Arg and the COM of the aromatic plane. | float | Å |




## 12. $\pi$ hydrogen bond

**pi_hbond**(trajectory, res_index_A, res_index_B, *frame=0, MIN_angle=120.0, MAX_distance_X_COM=5.5, MAX_distance_H_COM=3.0, MAX_distance_Hp_COM=1.2*)

### Description

Identify an hydrogen bond involving a $\pi$ area of an aromatic cycle.

It don't take in acount H of aromatic cycle, due to a false positive with aromatic-aromatic interaction, and the H of chaged residues (only the charged part), due to false positive with charge-aromatic interaction. 

Please note that protonated histidine are not taken in acount are not taken into account because they are involved in charge-aromatic interactions.

> [!IMPORTANT]
> The structure, or a trajectory must contain hydrogens.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| MIN_angle   | float | Minimum angle of X-H...COM. <br/> Default value: 120.0˚ | optional |
| MAX_distance_X_COM | float | Maximum distance between H donnor (X) and COM. <br/> Default value: 5.5 Å | optional |
| MAX_distance_H_COM | float | Maximum distance between H and COM <br/> Default value: 3.0 Å | optional |
| MAX_distance_Hp_COM | float | Maximum distance between H projected on aromatic plane and COM <br/> Default value: 1.2 Å | optional |


### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ), string (type), string (subtype) |  |
| .get_atoms    | Return a list of atoms index involved in H-bond.  | list of list: [[H_index, X_index],...] <br/> *X atom is the hydrogen donor.*|  |
| .get_angle    | List of angle betwee the N-COM and the normal of the aromatic ring. | list: [[angle],...] | degree |
| .get_distance | List of distance between COM of aromatic ring and N of the amino group. | list of list: [[distance_H_COM, distance_COM_Hprojected, distance_X_COM],...] | Å |

> [!TIP]
> In each list, the angle and distances have the same index as the list of atoms (H,X) making the $\pi$ H-bond.






## 13. $n \rightarrow \pi^*$

**n_pi**(trajectory, res_index_A, res_index_B, *frame=0, ref_distance=3.0, distance_tolerance=0.25, ref_angle=110 , angular_tolerance=10.0*)

### Description

Identify $n \rightarrow \pi^*$ between caronyl group (C=O) of two residues.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| ref_distance | float | Distance between O..C atoms. <br/> Default value: 3.0 Å | optional |
| distance_tolerance | float | The range of distance is ref_distance +/- N. <br/> Default value of N: 0.25 | optional |
| ref_angle | float | Angle between O..C=O atoms. <br/> Default value: 110.0 | optional |
| angular_tolerance | float | The range of angle is ref_angle +/- N. <br/> Default value of N: 10.0 | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ), string (regular / reciprocal) |  |
| .get_distance      | Return the two O...C distances. | floats: distance_O_res_A_C_res_B, distance_O_res_B_C_res_A | Å |
| .get_angle         | Return the two O...C=O angles.  | floats: angle_O_res_A_CO_res_B, angle_O_res_B_CO_res_A | degree |






## 14. Chalcogen bond & S/Se mediated H-bond

**sse_hydrogen_chalcogen_bond**(trajectory, res_index_A, res_index_B, *frame=0, MAX_distance_SX=3.6, MIN_angle_theta_chalcogen = 50.0, MAX_angle_dihedral_chalcogen=50.0, MIN_angle_phi_chalcogen=30.0, MAX_angle_phi_chalcogen=60.0, MIN_angle_csx_chalcogen=115.0, MAX_angle_csx_chalcogen=155.0, MAX_angle_csx_hbond=145.0*)   

### Description

Identify chalcogen bond or S/Se mediated H-bond. The analysis don't take in account the presence (or abscence) of hydrogens

> [!WARNING]
> The error `IndexError: list index out of range` Always happend when bonded atoms are missing. Ensure to have CONNECT section in your PDB, with the correct bonds.

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory  | mdtraj  | MDTraj trajectory.  | mandatory |
| res_index_A | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B | integer | Index of residue B in MDTraj topology. | mandatory |
| frame       | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| MAX_distance_SX              | integer | Cut-off distance between S/Se and the N/O to identify the interaction. <br/> Default value: 3.6 <br/> Unit: Angstrom | optional |
| MAX_angle_dihedral_chalcogen | integer | Maximum dihedral angle to identify chalcogen interaction. <br/> Default value: 50.0 <br/> Unit: degree | optional |
| MIN_angle_phi_chalcogen      | integer | Minimum phi angle to identify chalcogen interaction. <br/> Default value: 30.0 <br/> Unit: degree | optional |
| MAX_angle_phi_chalcogen      | integer | Maximum phi angle to identify chalcogen interaction. <br/> Default value: 60.0 <br/> Unit: degree | optional |
| MIN_angle_theta_chalcogen    | integer | Maximum theta angle to identify chalcogen interaction. <br/> Default value: 50.0 <br/> Unit: degree | optional |
| MIN_angle_csx_chalcogen      | integer | Minimum centroid-S/Se-N/O (CSX) angle to identify chalcogen interaction. <br/> Default value: 115.0 <br/> Unit: degree | optional |
| MAX_angle_csx_chalcogen      | integer | Maximum CSX angle to identify chalcogen interaction. <br/> Default value: 155.0 <br/> Unit: degree | optional |
| MAX_angle_csx_hbond          | integer | Maximum CSX angle to identify H-bond interaction. <br/> Default value: 145.0 <br/> Unit: degree | optional |

### Properties

| Property | Description | Return | Unit |
| -------- | --- | --- | --- |
| .check_interaction | Check if the given interaction type exisit.  | Boolean (True / False ) , string |  |
| .get_atoms      | Return list of atoms index involved in chalogen bond and another list of atoms index involved in H-bond. | List, List <br/> Format: [[index_S/Se,index_N/O],...] , [[index_S/Se,index_N/O],...] | MDTraj indices |
| .get_angle      | Return a list of angles for chalcogen bonds and another for H-bonds. | List, List <br/> Format: [[angle_theta, angle_phi, angle_CSX, angle_dihedral],...], [[angle_theta, angle_phi, angle_CSX, angle_dihedral],...] | Å |
| .get_distance   | Return a list of distance  for chalcogen bonds and another for H-bonds. | List, List <br/> Format: [distance,...] , [distance,...] | degree |





## 15. S/Se - Aromatic

**sse_aromatic**(trajectory, res_index_A, res_index_B, *frame=0, MAX_distance=5.5, MIN_pi_angle=60.0, MAX_quadrupole_angle=35.0*)

### Description

Identify the interaction between an sulfur or selenium atomand an aromatic resiue. The output give in which aromatic area the interaction take place: $\pi$, quadrupole or intermediate. 

### Arguments

| Argument | Format | Description | Requirement |
| -------- | --- | --- | --- |
| trajectory   | mdtraj  | MDTraj trajectory.  | mandatory |
| res_index_A  | integer | Index of residue A in MDTraj topology. | mandatory |
| res_index_B  | integer | Index of residue B in MDTraj topology. | mandatory |
| frame        | integer | Frame ID on which to perform the analysis. <br/> Default value: 0 | optional |
| MAX_distance | integer | Maximum distance between COM of aromatic ring and the S/Se. <br/> Default value: 5.5Å | optional |
| MIN_pi_angle | integer | Minimum angle defining the Pi area (the maximum is 90˚). <br/> Default value: 60.0˚ | optional |
| MAX_quadrupole_angle | integer | Maximum angle defining the quadrupole area (the maximum is 0˚). <br/> Default value: 35.0˚ | optional |

