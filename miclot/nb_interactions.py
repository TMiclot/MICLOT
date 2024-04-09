#!/usr/bin/env python3

"""
  ___      ___   __     ______   ___        ______  ___________  
 |"  \    /"  | |" \   /" _  "\ |"  |      /    " \("     _   ") 
  \   \  //   | ||  | (: ( \___)||  |     // ____  \)__/  \\__/  
  /\\  \/.    | |:  |  \/ \     |:  |    /  /    ) :)  \\_ /     
 |: \.        | |.  |  //  \ _   \  |___(: (____/ //   |.  |     
 |.  \    /:  | /\  |\(:   _) \ ( \_|:  \\        /    \:  |     
 |___|\__/|___|(__\_|_)\_______) \_______)\"_____/      \__|

 Molecular InteraCtion anaLysis tOolkiTs
 _______________________________________
 
 This module is provide function to compute the non-bonding interactions between two residues.
 
 The non-bonding interactions are:
    - C5 H-bond
    - C-bond
    - Hydrophobic & Repulsion hydrophobe/hydrophile
    - Charges clash & charges repulsion
    - Salt bridges
    - Hydrogen bond
    - van der Waals
    - 
"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"





#=====================================================
#===== Import modules
#=====================================================

import numpy as np
from skspatial.objects import Plane, Vector, Point
import mdtraj as md





#=====================================================
#===== Class for C5 hydrogen bonds
#=====================================================

class C5_hydrogen_bond:
    def __init__(self, trajectory, res_index, frame=0, angular_tolerance=10.0, MAX_distance=2.7):
        """
        INTERACTION TYPE    C5 hydrogen bond
        SUBTYPE(S)          No

        DESCRIPTION
            Intraresidue H-bond between carbonyl O and amino N atoms of the backbone.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index      Index of residue.
        
        OPTIONAL ARGUMENTS
            angular_tolerance     Set an absolute tolerance parameter N. See documentation concerning. 'numpy.isclose'.
                                  The range of angle is 140.0 +/- N.
                                  Default value of N: 10.0
            frame                 Frame ID on which to perform the analysis.
                                  Default value: 0
            MAX_distance          Maximum distance between O and H, in angstrom.
                                  Default value: 2.7
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res = res_index
        self.angular_tolerance = angular_tolerance
        self.MAX_distance = MAX_distance
        
        
        #===== Get atoms index for residue =====
        self.atom_res = self.top.select(f"resid {self.res}")
        self.atom_O_res = self.top.select(f"resid {self.res} and name O")[0]
        self.atom_H_res = self.top.select(f"resid {self.res} and name H HN")[0] # AMBER use the name H / CHARMM use the name HN
        
        
        #==== Compute O---H distance =====
        self.distance = md.compute_distances(traj, [[self.atom_O_res, self.atom_H_res]])[0][0] *10
        
        
        #===== Compute Psy and Phi angle =====
        self.all_psi_atoms, self.all_psi_angle = md.compute_psi(self.traj)
        self.all_psi_angle = np.rad2deg(self.all_psi_angle[0])
        #
        self.all_phi_atoms, self.all_phi_angle = md.compute_phi(self.traj)
        self.all_phi_angle = np.rad2deg(self.all_phi_angle[0])

        # Select value for the residue only
        for j_psi in self.all_psi_atoms:
            # select commom atoms between j_psi ant atoms in the residues
            common_elements_psi = np.intersect1d(j_psi,self.atom_res)
            
            for i_phi in self.all_phi_atoms:
                # select commom atoms between j_psi ant atoms in the residues
                common_elements_phi = np.intersect1d(i_phi,self.atom_res)
                
                # Check if the number of elements in the two arrays are greater than 2.
                # It mean: only if there are at least three atoms belonging to the residue
                if common_elements_phi.size > 2 and common_elements_psi.size > 2:
                    # index of atoms is the same as index of their angle
                    j_psi_index = np.where(self.all_psi_atoms == j_psi)[0][0]
                    i_phi_index = np.where(self.all_phi_atoms == i_phi)[0][0]
                    # get the phi and psi angles
                    self.angle_phi_res = self.all_phi_angle[i_phi_index]
                    self.angle_psi_res = self.all_psi_angle[j_psi_index]
    
    
    #===== Return parameters =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the C5 H-bond interaction exist for the residue.
        
        RETURN
            True     The interaction exist.
            False    The interaction don't exist.
        """
        # Take the absolute value of the angle to ensure negative and positive values are considered as the same (ex: -80 is 80)
        if self.distance <= self.MAX_distance \
        and np.isclose(abs(self.angle_phi_res), 140.0, atol=self.angular_tolerance) \
        and np.isclose(abs(self.angle_psi_res), 140.0, atol=self.angular_tolerance): #take absolute value of the angle for the checking
            return True
        else:
            return False
        
    @property
    def get_angle(self):
        """
        DESCRIPTION    Values of Phi and Psi angle of the residue.
        RETURN         angle_phi, angle_psi
        UNIT           degree
        """
        return self.angle_phi_res, self.angle_psi_res

    @property
    def get_distance(self):
        """
        DESCRIPTION    Distance between atoms O and H of the residue.
        RETURN         distance
        UNIT           Angstrom
        """
        return self.distance




      
#=====================================================
#===== Class for C-bond
#=====================================================

class C_bond:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0):
        """
        INTERACTION TYPE    C-bond
        SUBTYPE(S)          No

        DESCRIPTION
            Interaction between C(sp3) with N or O or S.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
           
       OPTIONAL ARGUMENTS
            frame    Frame ID on which to perform the analysis.
                     Default value: 0
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        
        
        
        #===== Create a dictionnaries of Csp3 and carbonyl O atoms in all amino acids =====
        self.Csp3 = {
            "ALA": "CA CB",
            "ARG": "CA CB CG CD",
            "ASN": "CA CB",
            "ASP": "CA CB",
            "CYS": "CA CB",
            "GLN": "CA CB CG",
            "GLU": "CA CB CG",
            "GLY": "CA",
            "HIS": "CA CB",
            "ILE": "CA CB CG1 CG2 CD",
            "LEU": "CA CB CG CD1 CD2",
            "LYS": "CA CB CG CD CE",
            "MET": "CA CB CG CE",
            "PHE": "CA CB",
            "PRO": "CA CB CG CD",
            "SER": "CA CB",
            "THR": "CA CB CG2",
            "TRP": "CA CB",
            "TYR": "CA CB",
            "VAL": "CA CB CG1 CG2",
            }
        
        self.O_carbonyl = {
            "ALA": "O",
            "ARG": "O",
            "ASN": "O OD1",
            "ASP": "O OD1",
            "CYS": "O",
            "GLN": "O OE1",
            "GLU": "O OE1",
            "GLY": "O",
            "HIS": "O",
            "ILE": "O",
            "LEU": "O",
            "LYS": "O",
            "MET": "O",
            "PHE": "O",
            "PRO": "O",
            "SER": "O",
            "THR": "O",
            "TRP": "O",
            "TYR": "O",
            "VAL": "O"
            }
        

        #===== Get atoms index for residue A =====
        self.res_A_name = self.top.residue(self.res_A).name
        self.atoms_Csp3_res_A = self.traj.topology.select(f"resid {self.res_A} and name {self.Csp3[self.res_A_name]}")
        self.atoms_O_res_A = self.traj.topology.select(f"resid {self.res_A} and name {self.O_carbonyl[self.res_A_name]}")
        
        #===== Get atoms index for residue B =====
        self.res_B_name = self.top.residue(self.res_B).name
        self.atoms_Csp3_res_B = self.traj.topology.select(f"resid {self.res_B} and name {self.Csp3[self.res_B_name]}")
        self.atoms_O_res_B = self.traj.topology.select(f"resid {self.res_B} and name {self.O_carbonyl[self.res_B_name]}")
        
        #===== Search all C-bonds between residues A and B =====
        # Create a list to store all C-bonds find
        # The list look like [[set_atoms_index_making_Cbond, list_theta, distance, energy], ...]
        self.list_Cbond = []
        
        # Search C-bond involved O in residue A and Csp3 in residue B
        for index_O in self.atoms_O_res_A:
            for index_C in self.atoms_Csp3_res_B:
                # Do the search
                interaction = self._search_Cbond(index_O, index_C)
                # If the variables are not empty
                if interaction != None:
                    self.list_Cbond.append(interaction)
                
        
        # Search C-bond involved O in residue B and Csp3 in residue A
        for index_O in self.atoms_O_res_B:
            for index_C in self.atoms_Csp3_res_A:
                # Do the search
                interaction = self._search_Cbond(index_O, index_C)
                # If the variables are not empty
                if interaction != None:
                    self.list_Cbond.append(interaction)



    #===== Function to search  =====
    # This is an internal function that should not be used by the user.
    # 
    def _search_Cbond(self, index_O, index_C):
        # Compute distance between O and C
        # *10 is used to convert nm to angstrom
        distance = md.compute_distances(self.traj, [[index_O,index_C]])[0][0] *10
        
        # Select distance between 2.5 and 3.6 A
        if 2.5 <= distance and distance <= 3.6:
            
            # Get the O and C atoms in the topology
            atom_O = self.traj.topology.atom(index_O)
            atom_C = self.traj.topology.atom(index_C)
            
            # Create sets to store the bonded atoms to O and C
            bonded_atoms_to_O = set()
            bonded_atoms_to_C = set()
            
            # Search bonded atoms to O and C
            for bond in self.traj.topology.bonds:
                # Check if bonded atom to O is C
                if atom_O in bond and bond[0].element.symbol in ['O','C'] and bond[1].element.symbol in ['O','C']:
                    bonded_atoms_to_O.add(bond[0])
                    bonded_atoms_to_O.add(bond[1])
                # Check if bonded atom to C is C, N or O
                elif atom_C in bond and bond[0].element.symbol in ['O','N','C'] and bond[1].element.symbol in ['O','N','C']:
                    bonded_atoms_to_C.add(bond[0])
                    bonded_atoms_to_C.add(bond[1])
            
            # Remove atom_O and atom_C from their bonded_atoms sets
            bonded_atoms_to_O.remove(atom_O)
            bonded_atoms_to_C.remove(atom_C)
            
            # Check possible theta_1 and theta_2 angles
            for y in bonded_atoms_to_O:
                # angle C...O=C (theta1)
                theta_1_atoms_indices = [[atom_C.index,atom_O.index,y.index]] #y is the C bonded to the O
                theta_1 = md.compute_angles(self.traj, theta_1_atoms_indices)[0][0]
                theta_1 = np.rad2deg(theta_1) #convert angle in radian to degree
                    
                for z in bonded_atoms_to_C:
                    # angle Z-C...O (theta2)
                    theta_2_atoms_indices = [[z.index,atom_C.index,atom_O.index]] # z is the atom bonded to the C
                    theta_2 = md.compute_angles(self.traj, theta_2_atoms_indices)[0][0]
                    theta_2 = np.rad2deg(theta_2) #convert angle in radian to degree
                    
                    # Compute interaction energy
                    energy = np.exp(np.log(22.2) - 22.2 * 10**6 * np.exp(-distance / 0.184)) - 22.2
                
                    # Check if geometric and energetic parameters are correct
                    if 160 <= theta_1 and theta_1 <=180 and 160 <= theta_2 and theta_2 <=180 and -22.2 <= energy and energy <=-2:
                        # create a set of atoms making the interaction
                        set_atoms_index_making_Cbond = {i for i in theta_1_atoms_indices[0] + theta_2_atoms_indices[0]}
                        # create a list containing all theta values
                        list_theta = [theta_1, theta_2]
                        
                        # return index of atom making c-bond, angles avlues, C-O distance and energy of the interaction 
                        return [set_atoms_index_making_Cbond, list_theta, distance, energy]
    
    
    
    #===== Return results & Properties=====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the C-bond interaction exist between tow residues.
        
        RETURN
            If the interaction exist between the 2 residues (it can be 1 or more interactactions):
                True
            
            If the interaction don't exist between the 2 residues:
                False
                
        OPTIONAL ARGUMENTS
            frame    Frame ID on which to perform the analysis.
                     Default value: 0
        """
        if len(self.list_Cbond) !=0:
            return True
        else:
            return False
        
    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return a list of atoms index involved in C-bond.
        RETURN         List of set: [{indices},{indices},...]
        """
        return [i[0] for i in self.list_Cbond]

    @property
    def get_angle(self):
        """
        DESCRIPTION
            Return a list of angles of the C-bond, it also return the corresponding C-bond atom indices.
            Theta 1: angle C...O=C
            Theta 2: angle Z-C...O
            
        RETURN
            List of list: [[[theta1,theta2],[indices]],...]
        """
        return [[i[1],i[0]] for i in self.list_Cbond]
        
    
    @property    
    def get_distance(self):
        """
        DESCRIPTION    Return a list of distance of the C-bond,
                       it also return the corresponding C-bond atom indices. 
        RETURN         List of list: [[[distance],[indices]],...]
        UNIT           Angstrom
        """
        return [[i[2],i[0]] for i in self.list_Cbond]        
        
    @property    
    def get_energy(self):
        """
        DESCRIPTION    Return a list of energy of the C-bond,
                       it also return the corresponding C-bond atom indices. 
        RETURN         List of list: [[[energy],[indices]],...]
        UNIT           kJ/mol
        """
        return [[i[3],i[0]] for i in self.list_Cbond] 




      
#=====================================================
#===== Class for Hydrophobic & Repulsion hydrophobe/hydrophile
#=====================================================

class hydrophobic:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_distance_COM=5.0, MAX_distance_CA=9.5, MIN_distance_CA=3.8):
        """
        INTERACTION TYPE    Hydrophobic or Hydrophobic/Hydrophilic repulsion
        SUBTYPE(S)          hydrophobic, clash

        DESCRIPTION
            Interaction between tow hydrophobic or hydrophilic residues.
            Distance between their center of mass of the side chain must be <= 5.0 A.
            Residues set as hydrophobic are: ALA, CYS, ILE, LEU, MET, PHE, TRP, TYR, VAL.
            Residues set as hydrophilic are: ARG, ASN, ASP, GLN, GLU, LYS
            
            * Compare to Pommié, C. et al., J. Mol. Recognit., 17, 17-32 (2004),
              TYR is considered as hydrophobic, instead of neutral.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
            
       OPTIONAL ARGUMENTS
            frame               Frame ID on which to perform the analysis.
                                Default value: 0
            MAX_distance_COM    Maximum distance between the COM of the side chain, in angstrom.
                                Default value: 5.0
            MAX_distance_CA     Maximum distance between the CA, in angstrom.
                                Default value: 9.5
            MIN_distance_CA     Minimum distance between the CA, in angstrom.
                                Default value: 3.8
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        self.MAX_distance_COM = MAX_distance_COM
        self.MAX_distance_CA = MAX_distance_CA
        self.MIN_distance_CA = MIN_distance_CA
        
        #===== Create a list of hydrophobic amino acids =====
        self.aa_hydrophobic = ['ALA','CYS','ILE','LEU','MET','PHE','TRP','TYR','VAL']
        self.aa_hydrophilic = ['ARG', 'ASN', 'ASP', 'GLN', 'GLU', 'LYS']
        self.list_aa_hydrophobic_hydrophilic = self.aa_hydrophobic + self.aa_hydrophilic
        
        #===== Get residue A and B names =====
        self.res_A_name = self.top.residue(self.res_A).name
        self.res_B_name = self.top.residue(self.res_B).name
        
        #===== Get CA atoms index and distance for residue A and B =====
        self.atom_CA_res_A = self.top.select(f"resid {self.res_A} and name CA")[0]
        self.atom_CA_res_B = self.top.select(f"resid {self.res_B} and name CA")[0]
        self.distance_CA = md.compute_distances(traj, [[self.atom_CA_res_A, self.atom_CA_res_B]])[0][0] *10 # *10 to convert nm to angstrom
        
        #===== Get center of mass (COM) of residue A and B & the distance between them =====
        # Compute COM of side chaines
        self.COM_res_A = md.compute_center_of_mass(self.traj, select=f"resid {self.res_A} and sidechain")[0]
        self.COM_res_B = md.compute_center_of_mass(self.traj, select=f"resid {self.res_B} and sidechain")[0]
        
        # Calculate distance beween COM (= Norm of the vector betweent he two COM)
        self.vector_COM_COM = Vector.from_points(self.COM_res_A, self.COM_res_B)
        self.distance_COM = self.vector_COM_COM.norm() *10 # *10 to convert nm to angstrom
    
    
    
    #===== Return parameters =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the hydrophobic attraction/clash exist between the two residues.
        
        RETURN
            True, True, hydrophobic      The interaction exist, and is: hydrophobic
            True, False, clash           The interaction exist, and is: hydrophobic/hydrophilic repulsion
            False, False, None           The interaction don't exist.
        """
        # check if the residues are hydrophobic or hydrophilic.
        if (self.res_A_name not in self.list_aa_hydrophobic_hydrophilic) or (self.res_B_name not in self.list_aa_hydrophobic_hydrophilic):
            return False, False, None
        
        elif self.distance_COM <= self.MAX_distance_COM and self.MIN_distance_CA <= self.distance_CA and self.distance_CA <= self.MAX_distance_CA \
             and self.res_A_name in self.aa_hydrophobic and self.res_B_name in self.aa_hydrophobic:
            return True, True, "hydrophobic"
        
        elif self.distance_COM <= self.MAX_distance_COM and self.MIN_distance_CA <= self.distance_CA and self.distance_CA <= self.MAX_distance_CA \
             and ((self.res_A_name in self.aa_hydrophobic and self.res_B_name in self.aa_hydrophilic) or (self.res_A_name in self.aa_hydrophilic and self.res_B_name in self.aa_hydrophobic)):
            return True, False, "clash"
        
        else:
            return False, False, None
    
    
    @property
    def get_distance(self):
        """
        DESCRIPTION    Distances between CA-CA atoms and side chaine COM-COM of the two residues.
        RETURN         distance_CA, distance_COM
        UNIT           Angstrom
        """
        return self.distance_CA, self.distance_COM




      
#=====================================================
#===== Class for charges clash & charges repulsion
#=====================================================

class charge_clash_repulsion:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_distance_CA=13.0, MIN_distance_charges=5.0):
        """
        INTERACTION TYPE    clash/repulstion
        SUBTYPE(S)          clash, repulsion

        DESCRIPTION
            Clash or repulstion between residue with same charge.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
            
       OPTIONAL ARGUMENTS
            frame                   Frame ID on which to perform the analysis.
                                    Default value: 0
            MAX_distance_CA         Maximum distance between CA, in angstrom
                                    Default value: 13.0
            MIN_distance_charges    Minimum distance to consider as clash, in angstrom.
                                    Default value: 5.0
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        self.MAX_distance_CA = MAX_distance_CA
        self.MIN_distance_charges = MIN_distance_charges
        
        #===== Create a lists of positive and negative amino acids =====
        self.aa_positif = ['LYS','ARG','HIP','HSP']
        self.aa_negatif = ['ASP','GLU']
        self.aa_charged = self.aa_positif + self.aa_negatif
        
        #===== Create a dictionnary of atoms representig charge position in charged AA =====
        self.charged_atoms = {
            "ARG": "CZ",
            "LYS": "NZ",
            "HIP": "ND1",
            "HSP": "ND1",
            "ASP": "CG",
            "GLU": "CD"
            }
        
        #===== Get atoms index & name for residue A =====
        self.atom_CA_res_A = self.top.select(f"resid {self.res_A} and name CA")[0]
        self.res_A_name = self.top.residue(self.res_A).name
        self.charged_atoms_res_A = self.top.select(f"resid {self.res_A} and name {self.charged_atoms[self.res_A_name]}")[0]
        
        #===== Get atoms index & name for residue B =====
        self.atom_CA_res_B = self.top.select(f"resid {self.res_B} and name CA")[0]
        self.res_B_name = self.top.residue(self.res_B).name
        self.charged_atoms_res_B = self.top.select(f"resid {self.res_B} and name {self.charged_atoms[self.res_B_name]}")[0]
        
        #==== Compute CA-CA and X-X distances =====
        self.distance_CA = md.compute_distances(traj, [[self.atom_CA_res_A, self.atom_CA_res_B]])[0][0] *10
        self.distance_charge = md.compute_distances(traj, [[self.charged_atoms_res_A, self.charged_atoms_res_B]])[0][0] *10
 


    #===== Return parameters =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the charge clash/repulsion interaction exist between the two residues.
        
        RETURN
            True, True, "clash"         The interaction exist, and is: clash
            True, False, "repulsion"    The interaction exist, and is: repulsion
            False, False, None          The interaction don't exist.
                                        Can happend if one or both residues are not charged, or have opposite charge.
        """
        # # check if at least one of the residues is not charged, or if they have opposite charges.
        if (self.res_A_name not in self.aa_charged) or (self.res_B_name not in self.aa_charged) \
        or (self.res_A_name in self.aa_positif and self.res_B_name in self.aa_negatif) \
        or (self.res_A_name in self.aa_negatif and self.res_B_name in self.aa_positif):
            return False, False, None
        else:
            # if CA distance is less, or equal, 13.0 A and charge-charge distance is less, or equal, 5.0 A.
            if self.distance_CA <= self.MAX_distance_CA and self.distance_charge <= self.MIN_distance_charges:
                return True, True, "clash"
            # if CA distance is less, or equal, 13.0 A and charge-charge distance is greater, or equal, 5.0 A.
            elif self.distance_CA <= self.MAX_distance_CA and self.MIN_distance_charges < self.distance_charge:
                return True, False, "repulsion"
            else:
                return False, False, None
    
    
    @property
    def get_distance(self):
        """
        DESCRIPTION    Distances between CA-CA and charge-charge atoms of the two residues.
        RETURN         distance_CA, distance_charge
        UNIT           Angstrom
        """
        return self.distance_CA, self.distance_charge





#=====================================================
#===== Class for salt bridges
#=====================================================

class salt_bridge:
    def __init__(self, trajectory, res_index_A, res_index_B, method="baker_hubbard", frame=0, MAX_distance_CA=13.0, MIN_distance_charges=4.0, distance_cutoff=3.0, angle_cutoff=120):
        """
        INTERACTION TYPE    Salt bridge
        SUBTYPE(S)          no

        DESCRIPTION
            Strong electrostatic interaction involving an H-bond and an
            ionic bond between amino acids with opposite charges.
            
            Electrostatic interaction is given by distances.
            
            H-bonds can be computed using 3 diffrent methods
            (see below, or take look at 'Hydrogen Bonding' in
            MDTraj documentation). For the method 'baker_hubbard' the
            following parameters are used:
                distance_cutoff=0.30    nm
                angle_cutoff=120        degree
            
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
         
         
       OPTIONAL ARGUMENTS
            frame              Frame ID on which to perform the analysis.
                               Have no effect with the method 'baker_hubbard',
                               because it implementation analyse the whole trajectory.
                               Default value: 0 
                               
            method             Use 'baker_hubbard' or 'kabsch_sander'
                               or 'wernet_nilsson' method to detect H-bond.
                               Default value: baker_hubbard
                               
            distance_cutoff         In the 'baker_hubbard' method, the distance cutoff of Donor-H…Acceptor contact, in angtrom. 
            angle_cutoff            In the 'baker_hubbard' method,the angle cutoff of the angle theta, in degrees. 
            MAX_distance_CA         Maximum distance between two CA atoms, in angtrom.
            MIN_distance_charges    Minimum distance between two charges, in angtrom.
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        self.method = method
        self.distance_cutoff = distance_cutoff/10
        self.angle_cutoff = angle_cutoff
        self.MAX_distance_CA = MAX_distance_CA
        self.MIN_distance_charges = MIN_distance_charges
        
        #===== Create a lists of positive and negative amino acids =====
        self.aa_positif = ['LYS','ARG','HIP','HSP']
        self.aa_negatif = ['ASP','GLU']
        self.aa_charged = self.aa_positif + self.aa_negatif
        
        #===== Create a dictionnary of atoms representig charge position in charged AA =====
        self.charged_atoms = {
            "ARG": "CZ",
            "LYS": "NZ",
            "HIP": "ND1", # Protonated HIS
            "HSP": "ND1", # Protonated HIS in CHARMM
            "ASP": "CG",
            "GLU": "CD"
            }
        
        #===== Create a dictionnary of atoms able to perform Hbond =====
        self.hbond_atoms = {
            "ARG": "HE NE NH1 HH11 HH12 NH2 HH21 HH22",
            "LYS": "NZ HZ1 HZ2 HZ3",
            "HIP": "ND1 NE2 HE2 HD1", # Protonated HIS
            "HSP": "ND1 NE2 HE2 HD1", # Protonated HIS in CHARMM
            "ASP": "OD1 OD2 HD2",
            "GLU": "OE1 OE2 HE2"
            }
        
        #===== Get atoms index & name for residue A =====
        self.atom_CA_res_A = self.top.select(f"resid {self.res_A} and name CA")[0]
        self.res_A_name = self.top.residue(self.res_A).name
        self.charged_atoms_res_A = self.top.select(f"resid {self.res_A} and name {self.charged_atoms[self.res_A_name]}")[0]
        
        #===== Get atoms index & name for residue B =====
        self.atom_CA_res_B = self.top.select(f"resid {self.res_B} and name CA")[0]
        self.res_B_name = self.top.residue(self.res_B).name
        self.charged_atoms_res_B = self.top.select(f"resid {self.res_B} and name {self.charged_atoms[self.res_B_name]}")[0]
        
        #==== Compute CA-CA and X-X distances =====
        # *10 to convert into Angstrom
        self.distance_CA = md.compute_distances(self.traj, [[self.atom_CA_res_A, self.atom_CA_res_B]])[0][0] *10
        self.distance_charge = md.compute_distances(self.traj, [[self.charged_atoms_res_A, self.charged_atoms_res_B]])[0][0] *10
        
        #==== Compute H-Bond =====
        # Create a new trajectory containing only the residue pair of atoms involve in H-bond of salt-bridge (other H-bond between the two residue are not considered)
        self.residue_pair = self.top.select(f"resid {self.res_A} {self.res_B} and name {self.hbond_atoms[self.res_A_name]} {self.hbond_atoms[self.res_B_name]}")
        self.traj_pair = self.traj.atom_slice(self.residue_pair)
        
        # Compute H-bond with selected method
        if self.method == "baker_hubbard":
            try:
                self.hbond = md.baker_hubbard(self.traj_pair, distance_cutoff=self.distance_cutoff, angle_cutoff=self.angle_cutoff, periodic=False) #done on the whole trajectory not a given frame
            except:
                self.hbond = np.array([])
        elif self.method == "kabsch_sander":
            try:
                self.hbond = md.kabsch_sander(self.traj_pair)
                self.hbond = self.hbond[0].toarray() #get only the matrix and convert it to numpy.array at the unique frame of 'traj_pair'
            except:
                self.hbond = np.array([])
        elif self.method == "wernet_nilsson":
            try:
                self.hbond = md.wernet_nilsson(self.traj_pair)
                self.hbond = self.hbond[0] #get only the numpy.array at selected at the unique frame of 'traj_pair'
            except:
                self.hbond = np.array([])

    
    #===== Return parameters =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the charge clash/repulsion interaction exist between the two residues.
        
        RETURN
            True, True      The interaction exist.
            True, False     Only charge attraction.
                            Can happend if the working file don't contain hydrogen atoms.
                            In this case, it's possible to considered that the interaction exists.
            False, True     Only H-bond between charged part of the residues.
            False, False    The interaction don't exist.
                            Can happend if one or both residues are not charged, or have same charge.
        """
        # check if at least one of the residues is not charged, or if both have the same charge.
        if (self.res_A_name not in self.aa_charged) or (self.res_B_name not in self.aa_charged) \
        or (self.res_A_name in self.aa_positif and self.res_B_name in self.aa_positif) \
        or (self.res_A_name in self.aa_negatif and self.res_B_name in self.aa_negatif):
            return False, False
        else:
            # if CA distance is less, or equal, 13 A and charge-charge distance is less or equal 4.0 A and the numpy.array is not empty.
            if self.distance_CA <= self.MAX_distance_CA and self.distance_charge <= self.MIN_distance_charges and np.count_nonzero(self.hbond) > 0:
                return True, True
            # if CA distance is less, or equal, 13 A and charge-charge distance is greater 4.0 A and the numpy.array is not empty.
            elif self.distance_CA <= self.MAX_distance_CA and self.MIN_distance_charges < self.distance_charge and np.count_nonzero(self.hbond) > 0:
                return False, True
            # if CA distance is less, or equal, 13 A and charge-charge distance is less or equal 4.0 A and the numpy.array is empty.
            elif self.distance_CA <= self.MAX_distance_CA and self.distance_charge <= self.MIN_distance_charges and np.count_nonzero(self.hbond) == 0:
                return True, False
            # other case
            else:
                return False, False    
    
    @property
    def get_distance(self):
        """
        DESCRIPTION    Distances between CA-CA and charge-charge atoms of the two residues.
        RETURN         distance_CA, distance_charge
        UNIT           Angstrom
        """
        return self.distance_CA, self.distance_charge

    
    @property
    def get_hbond(self):
        """
        DESCRIPTION    Return the atom indices forming H-bond in the pair.
                       For each method, the output is converted to an np.array.
        RETURN         hbond
        """ 
        
        if np.array_equal(self.hbond, np.array([])) or np.allclose(self.hbond, 0.):
            return np.array([])
        
        else:
            # Create a dictionnary of atom position and corresponding initial atom indeces
            dict_position_initial_indices = {" ".join(map(str,self.traj.xyz[0, atom_index])):atom_index for atom_index in self.residue_pair}

            # Create a dictionnary of atom indeces invlove in hbond and their position
            hbond_atom_indices = set(self.hbond.flatten())
            dict_hbond_indices_position = {atom_index:" ".join(map(str,self.traj_pair.xyz[0, atom_index])) for atom_index in hbond_atom_indices}

            # Create empty list to store all hbond
            list_hbond = []

            # Convert the 
            for i in self.hbond: #self.hbond is an np.array of np.array. Ex: [[0 2 3] [5 9 78] [65 2 7]]
                hbond = []
                for atom in i:
                    position = dict_hbond_indices_position[atom] #get the position of the atom in traj_pair
                    initial_index = dict_position_initial_indices[position] #use the position in traj_pair to get the initial atom index (positions are same in traj and traj_pair, only index change)
                    hbond.append(initial_index) # add the atom to the hbond interaction
                # Append list_hbond with the current hbond interaction, converted to np.array
                list_hbond.append(np.array(hbond)) 

            # Return the final np.array containig all hbonds
            return np.array(list_hbond)





#=====================================================
#===== Class for hydrogen bond
#=====================================================

class hydrogen_bond:
    def __init__(self, trajectory, res_index_A, res_index_B, method="baker_hubbard", distance_cutoff=3.0, angle_cutoff=150, frame=0):
        """
        INTERACTION TYPE    Hydrogen bond
        SUBTYPE(S)          regular, low-barrier, single-well,
                            intermediate_regular_low-barrier,
                            intermediate_low-barrier_single-well

        DESCRIPTION            
            H-bonds can be computed using 3 diffrent methods
            (see below, or take look at 'Hydrogen Bonding' in
            MDTraj documentation).
            This command don't compute H-bond involve in salt bridge,
            only salt_bridge do it.             
            
            For the method 'baker_hubbard' the following parameters are used:
                distance_cutoff=0.30    nm
                angle_cutoff=150        degree
            
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
            
       OPTIONAL ARGUMENTS
            frame     Frame ID on which to perform the analysis.
                      Have no effect with the method 'baker_hubbard',
                      because it implementation analyse the whole trajectory.
                      Default value: 0   
                      
            method    Use 'baker_hubbard' or 'kabsch_sander'
                      or 'wernet_nilsson' method to detect H-bond.
                      Default value: baker_hubbard
                      
            MAX_distance_CA         Maximum distance between two CA atoms, in angtrom.
            MIN_distance_charges    Minimum distance between two charges, in angtrom.
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        self.method = method
        self.distance_cutoff = distance_cutoff/10 #/10to convert angstrom to nm
        self.angle_cutoff = angle_cutoff
        
        #===== Create a lists of positive and negative amino acids =====
        self.aa_positif = ['LYS','ARG','HIP','HSP']
        self.aa_negatif = ['ASP','GLU']
        
        
        #===== Create a dictionnary of atoms that may be involved in a H-bond of a salt bridge =====
        self.hbond_atoms = {
            "ARG": "HE NE NH1 HH11 HH12 NH2 HH21 HH22",
            "LYS": "NZ HZ1 HZ2 HZ3",
            "HIP": "ND1 NE2 HE2 HD1", # Protonated HIS
            "HSP": "ND1 NE2 HE2 HD1", # Protonated HIS in CHARMM
            "ASP": "OD1 OD2 HD2",
            "GLU": "OE1 OE2 HE2"
            }
        
        
        #===== Get atoms index & name for the pair =====
        self.res_A_name = self.top.residue(self.res_A).name
        self.res_B_name = self.top.residue(self.res_B).name
        
        
        #===== identify Hbonds =====
        #if the residues are able to make a salt bridge, don't performd H-bond identification
        if (self.res_A_name in self.aa_positif and self.res_B_name in self.aa_negatif) \
        or (self.res_A_name in self.aa_negatif and self.res_B_name in self.aa_positif):
            self.hbond = np.array([]) 
            self.traj_pair = None
        
        # else perform H-bond identification
        else:
            # select the two residues forming the pair in the topology
            self.atoms_pair = self.top.select(f"resid {self.res_A} {self.res_B}")                    
                
            # Create a new trajectory containing only the tow residues
            self.traj_pair = self.traj.atom_slice(self.atoms_pair)

            # Compute all H-bond with selected method present in the system
            if self.method == "baker_hubbard":
                self.hbond = md.baker_hubbard(self.traj_pair, distance_cutoff=self.distance_cutoff, angle_cutoff=self.angle_cutoff, periodic=False) #done on the whole trajectory not a given frame
            elif self.method == "kabsch_sander":
                self.hbond = md.kabsch_sander(self.traj_pair)
                self.hbond = self.hbond[0].toarray() #get only the matrix and convert it to numpy.array at selected frame
            elif self.method == "wernet_nilsson":
                self.hbond = md.wernet_nilsson(self.traj_pair)
                self.hbond = self.hbond[0] #get only the numpy.array at selected frame
                
                
        #===== Convert output to np array =====
        if np.array_equal(self.hbond, np.array([])) or np.allclose(self.hbond, 0.):
            self.atom_indices = np.array([])
            
        else:
            # Create a dictionnary of atom position and corresponding initial atom indeces
            self.dict_position_initial_indices = {" ".join(map(str,self.traj.xyz[0, self.atom_index])):self.atom_index for self.atom_index in self.atoms_pair}

            # Create a dictionnary of atom indeces invlove in hbond and their position
            self.hbond_atom_indices = set(self.hbond.flatten())
            self.dict_hbond_indices_position = {self.atom_index:" ".join(map(str,self.traj_pair.xyz[0, self.atom_index])) for self.atom_index in self.hbond_atom_indices}

            # Create empty list to store all hbond
            self.list_bond = []

            # Previously we use 'traj.atom_slice' to create a traj containing only the pair.
            # But it change atom indices. So heres the correct atom indeces are assiged comparing the position.
            # If atom 0 have the same xzy position as atom 12 in the original traj, so the index 0 is replace by 12.
            for self.i in self.hbond: #self.hbond is an np.array of np.array. Ex: [[0 2 3] [5 9 78] [65 2 7]]
                self.bond = []
                for self.atom in self.i:
                    self.position = self.dict_hbond_indices_position[self.atom] #get the position of the atom in traj_pair
                    self.initial_index = self.dict_position_initial_indices[self.position] #use the position in traj_pair to get the initial atom index (positions are same in traj and traj_pair, only index change)
                    self.bond.append(self.initial_index) # add the atom to the hbond interaction
                # Append list_hbond with the current hbond interaction, converted to np.array
                self.list_bond.append(np.array(self.bond)) 

            # Return the final np.array containig all hbonds
            self.atom_indices = np.array(self.list_bond)
            

            
        #===== Compute distances and their energy profile (based on distance) =====
        # Create empty list to store all hbond and their properties
        self.list_Hbond = []
        # Initialise subtypes
        self.subtype = None
        
        for self.i in self.atom_indices:
            self.angle = np.rad2deg( md.compute_angles(self.traj, [self.i])[0][0] )
            self.distance_DA = md.compute_distances(self.traj, [[self.i[0], self.i[2]]])[0][0] *10
            if 2.7 <= self.distance_DA:
                self.subtype = "regular"
            elif 2.6 < self.distance_DA and self.distance_DA < 2.7:
                self.subtype = "intermediate_regular_low-barrier"
            elif 2.4 <= self.distance_DA and self.distance_DA <= 2.6:
                self.subtype = "low-barrier"
            elif 2.3 <= self.distance_DA and self.distance_DA < 2.4:
                self.subtype = "intermediate_low-barrier_single-well"
            elif self.distance_DA < 2.3:
                self.subtype = "single-well"
            
            self.list_Hbond.append([self.angle,self.distance_DA,self.subtype,self.i])
            
            
    #===== Return parameters =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the hydrogen bond interaction exist between the two residues.
        
        RETURN
            True      The interaction exist.
            False     The interaction don't exist.
                      Can happend if the working file don't contain hydrogen atoms,
                      or when it invole H-bond of a salt bridge.
        """
        # if the numpy.array is not empty.
        if np.count_nonzero(self.hbond) > 0:
            return True
        # if the numpy.array is empty.
        elif np.count_nonzero(self.hbond) == 0:
            return False
        # other case
        else:
            raise ValueError("MIN_contact_numbers must be greater or equal to 1.")    

    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return a list of atoms index involved in H-bond.
        RETURN         np.array : [[indices],[indices],...]
        """
        return self.atom_indices
        
    @property
    def get_angle(self):
        """
        DESCRIPTION    Return a list of angles (DHA) of all H-bonds,
                       it also return the corresponding atom indices.
        UNIT           Degree
        RETURN         List: [[angle,[indices]],...]
        """
        return [[i[0],i[3]] for i in self.list_Hbond]
        
    
    @property    
    def get_distance(self):
        """
        DESCRIPTION    Return a list of distance between donnor and acceptor,
                       it also return the corresponding C-bond atom indices. 
        RETURN         List: [[distance,[indices]],...]
        UNIT           Angstrom
        """
        return [[i[1],i[3]] for i in self.list_Hbond]

    
    @property    
    def get_subtype(self):
        """
        DESCRIPTION    Return the subtype of the H-bond. 
        RETURN         List: [[subtype,[indices]],...]
        UNIT           no
        """
        return [[i[2],i[3]] for i in self.list_Hbond]    
    
    



#=====================================================
#===== Class for van der Waals interaction
#=====================================================

class van_der_waals:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, set_hydrogen=True, distance_tolerance=0.5, MIN_contact_numbers=1):
        """
        INTERACTION TYPE    van der Waals
        SUBTYPE(S)          no

        DESCRIPTION
            Van der Waals interaction between residues.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS
            frame                 Frame ID on which to perform the analysis.
                                  Default value: 0
                                  
            distance_tolerance    The range of distance is radii_vdw_1 + radii_vdw_2 + N
                                  Range of values: 0.0 to 0.6 angstrom.
                                  Default value: 0.5
            
            set_hydrogen          Set if hydrogens are taken in acount for vdW interaction.
                                  Take a boolean: True or False.
                                  Default value: True
            
            MIN_contact_numbers   Min number of contacts tow residues must have to consider a vdW interaction.
                                  Default value: 1
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        self.set_hydrogen = set_hydrogen
        self.distance_tolerance = distance_tolerance
        self.MIN_contact_numbers = MIN_contact_numbers
        
        
        #===== Set if H are or not taken in account for vdw interaction =====
        # H is taken into account for vdW interaction
        if self.set_hydrogen == True:
            self.list_vdw_atoms = ["H","C","N","O","Se","S"]
        
        # H is NOT taken into account for vdW interaction
        elif self.set_hydrogen == False:
            self.list_vdw_atoms = ["C","N","O","Se","S"]
            
        
        #===== Create dictionarry of atoms radii =====
        # van der waals radii
        self.dict_radii_vdw = {
            "C":  1.70,
            "H":  1.10,
            "N":  1.55,
            "O":  1.52,
            "Se": 1.90, 
            "S":  1.80 
        }
        
        # covalent radii
        self.dict_radii_covalent = {
            "C":  0.75,
            "H":  0.32,
            "N":  0.71,
            "O":  0.64,
            "Se": 1.18, 
            "S":  1.04 
        }
        
        
        #===== Get residues informations =====
        self.res_A_topology = self.top.residue(self.res_A)
        self.res_B_topology = self.top.residue(self.res_B)
          
            
        #===== Compute VDW and interfaces between residues ====
        # VDW
        self.list_distance, self.list_contacts = self.__search_vdw()
        
        # interface
        self.interface = self.__get_interface()    
        
        
        
    #===== Function to search vdw =====
    # This is an internal function that should not be used by the user.
    # 
    def __search_vdw(self):        
        if self.distance_tolerance < 0 or self.distance_tolerance > 0.6:
            raise ValueError("Distance_tolerance must be between 0 and 0.6 anstrom.")
        else:
            #ceate list to set all distance and another with corresponding atom pairs
            list_distance = []
            list_contacts = []
            
            # loop over all atom of residue A
            for atom_A in self.res_A_topology.atoms:
                atom_A_symbol = atom_A.element.symbol
                
                # loop over all atom of residue B
                for atom_B in self.res_B_topology.atoms:
                    atom_B_symbol = atom_B.element.symbol
                    
                    # check if the symbol of the eleemtn of the two atoms are in the self.list_vdw_atoms
                    if atom_A_symbol in self.list_vdw_atoms and atom_B_symbol in self.list_vdw_atoms:
                        #compute covalent and vdw distance for the atom pair
                        distance_vdw = self.dict_radii_vdw[atom_A_symbol] + self.dict_radii_vdw[atom_B_symbol] + self.distance_tolerance
                        distance_covalent = self.dict_radii_covalent[atom_A_symbol] + self.dict_radii_covalent[atom_B_symbol]

                        # get distance between atoms. Use *10 to convert distance in nm to A
                        distance = md.compute_distances(self.traj, [[atom_A.index, atom_B.index]])[0][0] *10
                        
                        # check if the distance between atom_A and atom_B is grether thant their covalent distance
                        #     and lower or equal to their VDW distance + distance_tolerance
                        if distance_covalent < distance and distance <= distance_vdw:
                            list_distance.append(distance)
                            list_contacts.append([atom_A,atom_B])
    
        # return the result
        return list_distance, list_contacts
    
    
    
    #===== Function to compute interface between residues =====
    def __get_interface(self):
        #===== Computre interface between residues =====
        
        # H is taken into account for vdW interaction
        if self.set_hydrogen == True:
            sel_residue_pair = self.top.select(f"resid {self.res_A} {self.res_B}")
            sel_residue_A = self.top.select(f"resid {self.res_A}")
            sel_residue_B = self.top.select(f"resid {self.res_B}")
        elif self.set_hydrogen == False:
            sel_residue_pair = self.top.select(f"resid {self.res_A} {self.res_B} and not element H")
            sel_residue_A = self.top.select(f"resid {self.res_A} and not element H")
            sel_residue_B = self.top.select(f"resid {self.res_B} and not element H")
        
        # Isolate pair A-B from traj
        traj_pair = self.traj.atom_slice(sel_residue_pair)
        # Isolate residue A from traj
        traj_res_A = self.traj.atom_slice(sel_residue_A)
        # Isolate residue B from traj
        traj_res_B = self.traj.atom_slice(sel_residue_B)
        
        
        # Compute SASA of paire and residue
        #      probe_radius is in nm, here 0.001nm = 1pm
        #      n_sphere_points is the number of points representing the surface of each atom, here 1000
        SASA_pair  = md.shrake_rupley(traj_pair,  probe_radius=0.001, n_sphere_points=1000, mode='residue', change_radii=None, get_mapping=False)[0]
        SASA_res_A = md.shrake_rupley(traj_res_A, probe_radius=0.001, n_sphere_points=1000, mode='residue', change_radii=None, get_mapping=False)[0]
        SASA_res_B = md.shrake_rupley(traj_res_B, probe_radius=0.001, n_sphere_points=1000, mode='residue', change_radii=None, get_mapping=False)[0]
        
        # Compute interface
            # np.sum(self.SASA_pair) sum the SASA of the tow residue when forming pair
            # *100 to convert nm^2 to angstro^2
        interface = ( SASA_res_A[0] + SASA_res_B[0] - np.sum(SASA_pair) ) *100        
        
        # return result
        return interface
    
    
    #===== Return properties =====    
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the van der Waals interaction exist between the two residues.
        
        RETURN
            True      The interaction exist.
            False     The interaction don't exist.
        """
        if self.MIN_contact_numbers <=0:
            raise ValueError("MIN_contact_numbers must be greater or equal to 1.")
        else:
            # if the list of contact contain a number equal or gereter of element, return true 
            if self.MIN_contact_numbers <= len(self.list_contacts): 
                return True
            else:
                return False
     
        
    @property
    def get_distance(self):
        """
        DESCRIPTION    Return the list distances between atoms pairs making vdW interaction and the list of their index.
        RETURN         list_distance, list_contacts
        UNIT           Angstrom
        """
        return self.list_distance, self.list_contacts
    
        
    @property
    def get_number_contacts(self):
        """
        DESCRIPTION    Return the number of atom-atom vdW contact between the two residues.
        RETURN         size of the list_contacts
        """
        return len(self.list_contacts)

    
    @property
    def get_interface(self):
        """
        DESCRIPTION    Return the interface contact between the two residues.
        RETURN         interface of residue A, interface of residue B
                       None is return when the vdw interaction don't exist
        UNIT           angstrom^2
        """
        #check if the vdw interaction exist between the tow residue, if es return the interface else return None
        if 0 < len(self.list_contacts):
            return self.interface
        else:
            return None
        




#=====================================================
#===== Class for 
#=====================================================

