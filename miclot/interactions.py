#!/usr/bin/env python3

"""
This module is provide function to compute the non-bonding interactions between two residues.

"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"


__all__ = ['C5_hydrogen_bond', 'C_bond', 'hydrophobic', 'charge_clash_repulsion', 'salt_bridge', \
           'hydrogen_bond', 'van_der_waals', 'amino_pi', 'charge_aromatic', 'aromatic_aromatic', \
           'arg_involved', 'pi_hbond', 'n_pi', 'sse_hydrogen_chalcogen_bond', 'sse_aromatic', \
           'identify_all_interaction_pair', 'interaction_table_whole_system']



#=====================================================
#===== Import modules
#=====================================================

import os
import glob
from collections import Counter

import numpy as np
from skspatial.objects import Plane, Vector, Points
import mdtraj as md
import pandas as pd

from multiprocessing import Pool

from itertools import combinations
from tqdm import tqdm

from .coulomb_lj import omm_coulomb_lj


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
        self.distance = md.compute_distances(self.traj, [[self.atom_O_res, self.atom_H_res]])[0][0] *10
        
        
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
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         [[atom_O_res, atom_H_res]]
        """
        return [[self.atom_O_res, self.atom_H_res]]


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
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_distance=3.6, MIN_distance=2.5, \
                 MAX_angle_COC=180.0, MIN_angle_COC=160.0, MAX_angle_ZOC=180.0, MIN_angle_ZOC=160.0, \
                 MAX_energy=-2.0, MIN_energy=-22.2):
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
                     
            MAX_distance   Maximum C-O distance
                           Default value: 3.6 angstrom
            MIN_distance   Mimimum C-O distance
                           Default value: 2.5 angstrom
            
            MAX_angle_COC    Maximum C..O=C angle
                             Default value: 180.0 degree
            MIN_angle_COC    Minimum C..O=C angle
                             Default value: 160.0 degree
            
            MAX_angle_ZOC    Maximum Z-O...C angle
                             Default value: 180.0 degree
            MIN_angle_ZOC    Minimum Z-O...C angle
                             Default value: 160.0 degree
            
            MAX_energy    Maximum energy of the interaction
                          Default value: -2.0 kJ/mol
            MIN_energy    Minimum energy of the interaction
                          Default value: -22.0 kJ/mol
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        self.MAX_distance  = MAX_distance
        self.MIN_distance  = MIN_distance
        self.MAX_angle_COC = MAX_angle_COC
        self.MIN_angle_COC = MIN_angle_COC
        self.MAX_angle_ZOC = MAX_angle_ZOC
        self.MIN_angle_ZOC = MIN_angle_ZOC
        self.MAX_energy = MAX_energy
        self.MIN_energy = MIN_energy
        
        
        
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
        self.res_A_name = self.traj.topology.residue(self.res_A).name
        self.atoms_Csp3_res_A = self.traj.topology.select(f"resid {self.res_A} and name {self.Csp3[self.res_A_name]}")
        self.atoms_O_res_A = self.traj.topology.select(f"resid {self.res_A} and name {self.O_carbonyl[self.res_A_name]}")
        
        #===== Get atoms index for residue B =====
        self.res_B_name = self.traj.topology.residue(self.res_B).name
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
        
        # Select distance between MIN_ MAX_ distances
        if self.MIN_distance <= distance <= self.MAX_distance:
            
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
                    if self.MIN_angle_COC <= theta_1 <= self.MAX_angle_COC and self.MIN_angle_ZOC <= theta_2 <= self.MAX_angle_ZOC and self.MIN_energy <= energy <= self.MAX_energy:
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
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_distance_COM=5.0, MAX_distance_CA=9.5):
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
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        self.MAX_distance_COM = MAX_distance_COM
        self.MAX_distance_CA = MAX_distance_CA
        
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
        self.distance_CA = md.compute_distances(self.traj, [[self.atom_CA_res_A, self.atom_CA_res_B]])[0][0] *10 # *10 to convert nm to angstrom
        
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
            True, False, clash           The interaction exist, and is: hydrophobic/hydrophilic clash
            True, False, repulsion       The interaction exist, and is: hydrophobic/hydrophilic repulsion
            False, False, None           The interaction don't exist.
        """
        # check if the residues are hydrophobic or hydrophilic.
        if (self.res_A_name not in self.list_aa_hydrophobic_hydrophilic) or (self.res_B_name not in self.list_aa_hydrophobic_hydrophilic):
            return False, False, None
        
        elif self.distance_COM <= self.MAX_distance_COM and self.distance_CA <= self.MAX_distance_CA \
             and self.res_A_name in self.aa_hydrophobic and self.res_B_name in self.aa_hydrophobic:
            return True, True, "hydrophobic"
        
        elif self.distance_COM <= self.MAX_distance_COM and self.distance_CA <= self.MAX_distance_CA \
             and ((self.res_A_name in self.aa_hydrophobic and self.res_B_name in self.aa_hydrophilic) or (self.res_A_name in self.aa_hydrophilic and self.res_B_name in self.aa_hydrophobic)):
            return True, False, "clash"
        
        # if the COM-COM distance is more than the MAX_distance_COM, interactioinis considered as repulsion.
        elif self.MAX_distance_COM <= self.distance_COM and self.distance_CA <= self.MAX_distance_CA \
             and ((self.res_A_name in self.aa_hydrophobic and self.res_B_name in self.aa_hydrophilic) or (self.res_A_name in self.aa_hydrophilic and self.res_B_name in self.aa_hydrophobic)):
            return True, False, "repulsion"
        
        else:
            return False, False, None
    

    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         [atoms_interaction]
        """
        atoms_interaction = []
        atoms_interaction.extend(list(self.top.select(f"resid {self.res_A} and sidechain")))
        atoms_interaction.extend(list(self.top.select(f"resid {self.res_B} and sidechain")))
        return [atoms_interaction]


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
        self.distance_CA = md.compute_distances(self.traj, [[self.atom_CA_res_A, self.atom_CA_res_B]])[0][0] *10
        self.distance_charge = md.compute_distances(self.traj, [[self.charged_atoms_res_A, self.charged_atoms_res_B]])[0][0] *10
 


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
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         [atoms_interaction]
        """
        atoms_interaction = []
        atoms_interaction.extend([self.charged_atoms_res_A])
        atoms_interaction.extend([self.charged_atoms_res_B])
        return [atoms_interaction]


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
            "HIP": "ND1 NE2 HE2 HD1 HD2", # Protonated HIS
            "HSP": "ND1 NE2 HE2 HD1 HD2", # Protonated HIS in CHARMM
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
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         [atoms_interaction]
        """
        atoms_interaction = []
        atoms_interaction.extend([self.charged_atoms_res_A])
        atoms_interaction.extend([self.charged_atoms_res_B])
        atoms_interaction.extend(list(self.residue_pair)) # part of residue pair containing only atom performing Hbond in salt bridge
        return [atoms_interaction]
        
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
            "HIP": "ND1 NE2 HE2 HD1 HD2", # Protonated HIS
            "HSP": "ND1 NE2 HE2 HD1 HD2", # Protonated HIS in CHARMM
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
        SUBTYPE(S)          charged-charged, charged-polar, charged-apolar,
                            polar-polar, polar-apolar, apolar-apolar

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
        
        # Get names
        self.res_A_name = self.top.residue(self.res_A).name
        self.res_B_name = self.top.residue(self.res_B).name  
          
            
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
            # *100 to convert nm^2 to angstrom^2
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
            True, 'subtype'    The interaction exist.
            False, False       The interaction don't exist.
        """
        # List of AA properties
        polar = ['CYS', 'HIS', 'HIE', 'HID', 'HSE', 'HSD', 'ASN', 'GLN', 'SER', 'THR', 'TRP']
        apolar = ['ALA', 'PHE', 'GLY', 'ILE', 'LEU', 'VAL', 'MET', 'PRO', 'TYR']
        charged = ['GLU', 'ASP', 'LYS', 'ARG', 'HIP', 'HSP', 'CYM', 'TYM']
        
        if self.res_A_name in polar and self.res_B_name in polar:
            subtype = "polar-polar"
        elif self.res_A_name in apolar and self.res_B_name in apolar:
            subtype = "apolar-apolar"
        elif self.res_A_name in charged and self.res_B_name in charged:
            subtype = "charged-charged"
        elif (self.res_A_name in polar and self.res_B_name in apolar) or (self.res_B_name in polar and self.res_A_name in apolar):
            subtype = "polar-apolar"
        elif (self.res_A_name in charged and self.res_B_name in apolar) or (self.res_B_name in charged and self.res_A_name in apolar):
            subtype = "charged-apolar"
        elif (self.res_A_name in charged and self.res_B_name in polar) or (self.res_B_name in charged and self.res_A_name in polar):
            subtype = "charged-polar"
        
        
        # Check geomatric parameters
        if self.MIN_contact_numbers <=0:
            raise ValueError("MIN_contact_numbers must be greater or equal to 1.")
        else:
            # if the list of contact contain a number equal or gereter of element, return true 
            if self.MIN_contact_numbers <= len(self.list_contacts): 
                return True, subtype
            else:
                return False, False
     
        
    @property
    def get_distance(self):
        """
        DESCRIPTION    Return the list distances between atoms pairs making vdW interaction and the list of their index.
        RETURN         list_distance, list_contacts
        UNIT           Angstrom
        """
        return self.list_distance
    
    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return the list of atom pairs performing vdw contact
        RETURN         list of atom pair list as [[indices], [indices], ...]
        """
        return [[y.index for y in x] for x in self.list_contacts]
    
        
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
#===== Class for Amino-Pi
#=====================================================
class amino_pi:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_distance=5.5, angular_tolerance=30.0):
        """
        INTERACTION TYPE    amino-Pi
        SUBTYPE(S)          no

        DESCRIPTION
            Interaction the amino group (-NH2) of ASN, or GLN, and an aromatic rings.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS
        Set an absolute tolerance parameter N. See documentation concerning. 'numpy.isclose'.
        
            MAX_distance          Maximum distance between COM of aromatic ring and the N of the amino group
                                  Default value: 5.5

            angular_tolerance     The range of angle is 90.0 +/- N.
                                  Default value of N: 30.0

            frame                 Frame ID on which to perform the analysis.
                                  Default value: 0: 0
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.MAX_distance = MAX_distance
        self.angular_tolerance = angular_tolerance
        self.res_index_A = res_index_A
        self.res_index_B = res_index_B
        
        
        #===== Create dictionaries to  =====
        self.dict_aromatic_ring = {"TYR": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "TRP": 'CD2 CE2 CE3 CZ2 CZ3 CH2',
                                   "PHE": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "HIS": 'CG ND1 CD2 CE1 NE2',
                                   "HID": 'CG ND1 CD2 CE1 NE2',
                                   "HIE": 'CG ND1 CD2 CE1 NE2',
                                   "HSD": 'CG ND1 CD2 CE1 NE2',
                                   "HSE": 'CG ND1 CD2 CE1 NE2',
                                  }
        
        
        self.dict_plane = {"TYR": ['CG',  'CE1',  'CE2'],
                           "TRP": ['CD2', 'CZ2',  'CZ3'],
                           "PHE": ['CG',  'CE1',  'CE2'],
                           "HIS": ['CG',  'CE1',  'NE2'],
                           "HID": ['CG',  'CE1',  'NE2'],
                           "HIE": ['CG',  'CE1',  'NE2'],
                           "HSE": ['CG',  'CE1',  'NE2'],
                           "HSD": ['CG',  'CE1',  'NE2'],
                          }

        
        
        #===== Identify residues and get their names =====
        # Identify the aromatic and the amino AA
        if self.top.residue(self.res_index_A).name in self.dict_aromatic_ring and self.top.residue(self.res_index_B).name in ["ASN", "GLN"]:
            self.res_amino = self.res_index_B
            self.res_aromatic = self.res_index_A
        elif self.top.residue(self.res_index_B).name in self.dict_aromatic_ring and self.top.residue(self.res_index_A).name in ["ASN", "GLN"]:
            self.res_amino = self.res_index_A
            self.res_aromatic = self.res_index_B
        else:
            raise ValueError('Residues are not TYR TRP PHE HIS / ASN GLN')
            
        # Get residues names
        self.res_amino_name = self.top.residue(self.res_amino).name
        self.res_aromatic_name = self.top.residue(self.res_aromatic).name

        
        #===== Distance between COM of aromatic and N =====
        # Compute COM of aromatic ring
        self.COM_aromatic = md.compute_center_of_mass(self.traj, select=f"resid {self.res_aromatic} and name {self.dict_aromatic_ring[self.res_aromatic_name]}")[0]
        
        # Get position of the N
        self.atom_N_index = self.top.select(f"resid {self.res_amino} and name ND2 NE2")[0]
        self.atom_N_position = self.traj.xyz[0][self.atom_N_index]
        
        # create a vector between aromatic COM and N & measure it's norm
        self.vector_COM_N = Vector.from_points(self.COM_aromatic, self.atom_N_position)
        self.distance = self.vector_COM_N.norm() *10 # *10 to convert nm to angstrom
        
        
        #===== Create planes and get it's normal vector =====
        # Get position of all atom in aromatic plane
        self.list_aromatic_atom_position = []
        #
        for self.atom in self.dict_plane[self.res_aromatic_name]:
            self.atom_index = self.top.residue(self.res_aromatic).atom(self.atom).index
            self.atom_position = self.traj.xyz[0][self.atom_index]
            self.list_aromatic_atom_position.append(self.atom_position)
        
        # create aromatic plane
        self.plane_aromatic = Plane.from_points(self.list_aromatic_atom_position[0], self.list_aromatic_atom_position[1], self.list_aromatic_atom_position[2])
        
        # Get normal vector od the aromatic plane
        self.vector_normal_plane_aromatic = self.plane_aromatic.normal
        
         
        # ===== Measure angle between the normal vector of aromatic plane and the vector COM --> N =====
        # calculate the angle and convert it to degree
        self.angle = np.rad2deg( self.vector_normal_plane_aromatic.angle_between(self.vector_COM_N) )
        
        # Convert the angle to a range of values from 0 to 90 degrees
        if self.angle <= 90:
            None # the angle is already in the range 0 to 90 degrees
        elif self.angle <= 180:
            self.angle = 180 - self.angle
        elif self.angle <= 270:
            self.angle = self.angle - 180
        else:
            self.angle = 360 - self.angle
        
        # The angle is before is between the normal vector of the plan and the vector COM --> N
        # But the normal vector is perpendicular to the plan and geometric criteria are for the angle between the plan and the vector COM --> N.
        # So angle is converted to correspond to be the one betwen the plan and the vector COM --> N.
        self.angle = 90 - self.angle
    
    
    #===== Return results =====
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
        if self.distance <= self.MAX_distance and np.isclose(self.angle, 90.0, atol=self.angular_tolerance):
            return True
        else:
            return False

    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         [atoms_interaction]
        """
        atoms_interaction = list(self.top.select(f"resid {self.res_aromatic} and name {self.dict_aromatic_ring[self.res_aromatic_name]}"))
        atoms_interaction.append(self.atom_N_index)
        return [atoms_interaction]

    @property
    def get_angle(self):
        """
        DESCRIPTION    Angle betwee the N-COM and the normal of the amino plane
        RETURN         angle_planes
        UNIT           degree
        """
        return abs(self.angle)

    @property
    def get_distance(self):
        """
        DESCRIPTION    Distance between COM of aromatic ring and N of the amino group
        RETURN         distance
        UNIT           Angstrom
        """
        return self.distance
    
    







#=====================================================
#===== Class for charge-aromatic interactions
#=====================================================
class charge_aromatic:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_distance=5.5, MIN_pi_angle=60.0, MAX_quadrupole_angle=35.0):
        """
        INTERACTION TYPE    charge with aromatic ring
        SUBTYPE(S)          cation-pi, anion-pi / cation-intermediate, anion-intermediate / cation-quadrupole, anion-quadrupole

        DESCRIPTION
            Interaction with the charged amino acid and an aromatic rings.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS        
            MAX_distance          Maximum distance between COM of aromatic ring and the charge
                                  Default value: 5.5
                                  
            MIN_pi_angle          Minumum angle to set Pi area. (The max angle is 90)
                                  Default value: 60.0
                                  
            MAX_quadrupole_angle   Maximum angle to set quadrupole area. (The min angle is 0)
                                   Default value: 35.0

            frame                 Frame ID on which to perform the analysis.
                                  Default value: 0
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.MAX_distance = MAX_distance
        self.MIN_pi_angle = MIN_pi_angle
        self.MAX_quadrupole_angle = MAX_quadrupole_angle
        self.res_index_A = res_index_A
        self.res_index_B = res_index_B

        
        #===== Create dictionaries of aromatic ring and corresponding plans =====
        self.dict_aromatic_ring = {"TYR": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "TRP": 'CD2 CE2 CE3 CZ2 CZ3 CH2',
                                   "PHE": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "HIS": 'CG ND1 CD2 CE1 NE2',
                                   "HID": 'CG ND1 CD2 CE1 NE2',
                                   "HIE": 'CG ND1 CD2 CE1 NE2',
                                   "HSD": 'CG ND1 CD2 CE1 NE2',
                                   "HSE": 'CG ND1 CD2 CE1 NE2',
                                  }
           
        self.dict_plane = {"TYR": ['CG',  'CE1',  'CE2'],
                           "TRP": ['CD2', 'CZ2',  'CZ3'],
                           "PHE": ['CG',  'CE1',  'CE2'],
                           "HIS": ['CG',  'CE1',  'NE2'],
                           "HID": ['CG',  'CE1',  'NE2'],
                           "HIE": ['CG',  'CE1',  'NE2'],
                           "HSE": ['CG',  'CE1',  'NE2'],
                           "HSD": ['CG',  'CE1',  'NE2'],
                          }

        
        #===== Create a dictionnary of atoms representig charge position in charged AA =====
        self.dict_charged_atoms = {
            "ARG": "CZ",
            "LYS": "NZ",
            "HIP": "ND1", # Protonated HIS
            "HSP": "ND1", # Protonated HIS in CHARMM
            "ASP": "CG",
            "GLU": "CD"
            }
        
        
        #===== Identify residues and get their names =====
        # Identify the aromatic and the amino AA
        if self.top.residue(res_index_A).name in self.dict_aromatic_ring and self.top.residue(res_index_B).name in self.dict_charged_atoms:
            self.res_charge = res_index_B
            self.res_aromatic = res_index_A
        elif self.top.residue(res_index_B).name in self.dict_aromatic_ring and self.top.residue(res_index_A).name in self.dict_charged_atoms:
            self.res_charge = res_index_A
            self.res_aromatic = res_index_B
        else:
            raise ValueError('Residues are not TYR TRP PHE HIS / ARG LYS ASP GLU HIP(HSP)')
            
        # Get residues names
        self.res_charge_name = self.top.residue(self.res_charge).name
        self.res_aromatic_name = self.top.residue(self.res_aromatic).name
    
        
        #===== Distance between COM of aromatic and charge =====
        # Compute COM of aromatic ring
        self.COM_aromatic = md.compute_center_of_mass(self.traj, select=f"resid {self.res_aromatic} and name {self.dict_aromatic_ring[self.res_aromatic_name]}")[0]
        
        # Get position of the charge
        self.atom_charge_index = self.top.select(f"resid {self.res_charge} and name {self.dict_charged_atoms[self.res_charge_name]}")[0]
        self.atom_charge_position = self.traj.xyz[0][self.atom_charge_index]
        
        # create a vector between aromatic COM and charge & measure it's norm
        self.vector_COM_charge = Vector.from_points(self.COM_aromatic, self.atom_charge_position)
        self.distance = self.vector_COM_charge.norm() *10 # *10 to convert nm to angstrom
        
        
        #===== Create planes and get it's normal vector =====
        # Get position of all atom in aromatic plane
        self.list_aromatic_atom_position = []
        #
        for self.atom in self.dict_plane[self.res_aromatic_name]:
            self.atom_index = self.top.residue(self.res_aromatic).atom(self.atom).index
            self.atom_position = self.traj.xyz[0][self.atom_index]
            self.list_aromatic_atom_position.append(self.atom_position)
        
        # create aromatic plane
        self.plane_aromatic = Plane.from_points(self.list_aromatic_atom_position[0], self.list_aromatic_atom_position[1], self.list_aromatic_atom_position[2])
        
        # Get normal vector od the aromatic plane
        self.vector_normal_plane_aromatic = self.plane_aromatic.normal
        
         
        # ===== Measure angle between the normal vector of aromatic plane and the vector COM --> charge =====
        # calculate the angle and convert it to degree
        self.angle = np.rad2deg( self.vector_normal_plane_aromatic.angle_between(self.vector_COM_charge) )
        
        # Convert the angle to a range of values from 0 to 90 degrees
        if self.angle <= 90:
            None # the angle is already in the range 0 to 90 degrees
        elif self.angle <= 180:
            self.angle = 180 - self.angle
        elif self.angle <= 270:
            self.angle = self.angle - 180
        else:
            self.angle = 360 - self.angle
        
        # The angle is before is between the normal vector of the plan and the vector COM --> charge
        # But the normal vector is perpendicular to the plan and geometric criteria are for the angle between the plan and the vector COM --> charge.
        #
        #    Pn                Pn is the normal vector to Plane.
        #    |   v             We calculate the angle between Pn and v.
        #    |  /              We need the angle between V and Plane. 
        #    | /
        # ___|/_____ Plane
        #
        # So angle is converted to correspond to be the one betwen the plan and the vector COM --> charge.
        self.angle = 90 - self.angle
        
        
        
    #===== Return results =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if a charged residue interact with an aromatic ring.
        
        RETURN  
            When True, 'charge' is replaced by cation or anion.
            True, 'charge'-pi              The interaction exist.
            True, 'charge'-intermediate    The interaction exist.
            True, 'charge'-quadrupole      The interaction exist.
            
            False, False                   The interaction don't exist.
        """
        # get charge type
        if self.res_charge_name in ["ARG","LYS","HIP","HSP"]:
            charge_type = "cation"
        elif self.res_charge_name in ["GLU","ASP"]:
            charge_type = "anion"
        
        # Take the absolute value of the angle to ensure negative and positive values are considered as the same (ex: -80 is 80)
        if self.distance <= self.MAX_distance and self.MIN_pi_angle <= self.angle:
            return True, f"{charge_type}-pi"
        elif self.distance <= self.MAX_distance and self.MAX_quadrupole_angle < self.angle < self.MIN_pi_angle:
            return True, f"{charge_type}-intermediate"
        elif self.distance <= self.MAX_distance and self.angle <= self.MAX_quadrupole_angle:
            return True, f"{charge_type}-quadrupole"
        else:
            return False, False
    
    
    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         [atoms_interaction]
        """
        atoms_interaction = list(self.top.select(f"resid {self.res_aromatic} and name {self.dict_aromatic_ring[self.res_aromatic_name]}"))
        atoms_interaction.append(self.atom_charge_index)
        return [atoms_interaction]


    @property
    def get_angle(self):
        """
        DESCRIPTION    Angle betwee the charge-COM and the normal of the aromatic ring
        RETURN         angle_planes
        UNIT           degree
        """
        return abs(self.angle)


    @property
    def get_distance(self):
        """
        DESCRIPTION    Distance between COM of aromatic ring and N of the amino group
        RETURN         distance
        UNIT           Angstrom
        """
        return self.distance
    
    





#=====================================================
#===== Class for aromatic-aromatic interactions =====
#=====================================================    
class aromatic_aromatic:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_angle_planarity=30.0, MAX_distance_COM=5.5, MIN_distance_offset=1.6, MAX_distance_offset=2.0, \
                 MIN_pi_angle=60.0, MAX_quadrupole_angle=35.0, MAX_angle_Tshaped=5.0):
        """
        INTERACTION TYPE    Stacking of two aromatic residues.
        SUBTYPE(S)          parallel, offset, t-shaped, y-shaped, coplanar, intermediate

        DESCRIPTION
            Stacking of two aromatic rings.
            Please note that protonated histidine are not taken in acount are not taken into account because they can make charge-aromatic interactions.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS

            frame    Frame ID on which to perform the analysis.
                     Default value: 0
                     
            MAX_distance_COM        Distance between the COMs of each residue
                                    Default value: 5.5 Å
            MIN_distance_offset     Minimum distance between the COM of the first residue and the COM of the second residue projected on the plane of the first residue.
                                    Default value: 1.6 Å
            MAX_distance_offset     Maximum distance between the COM of the first residue and the COM of the second residue projected on the plane of the first residue.
                                    Default value: 2.0 Å
            MIN_pi_angle            Minimum angle defining the Pi area (the maximum is 90˚).
                                    Default value: 60.0˚
            MAX_quadrupole_angle    Maximum angle defining the quadrupole area (the maximum is 0˚).
                                    Default value: 35.0˚
            MAX_angle_planarity     Maximum angle defining the planarity between the two planes (the maximum is 0˚).
                                    Default value: 30.0˚
            MAX_angle_Tshaped       Maximum angle between the normal vector of an aromatic plane of residue 1 and the vector COM-->C of the ring of the residue 2
                                    (the maximum is 0˚).
                                    Default value: 5.0˚
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_index_A = res_index_A
        self.res_index_B = res_index_B
        self.MAX_distance_COM = MAX_distance_COM # distance between the COM-COM of the two residues
        self.MIN_distance_offset = MIN_distance_offset
        self.MAX_distance_offset = MAX_distance_offset
        self.MIN_pi_angle = MIN_pi_angle
        self.MAX_quadrupole_angle = MAX_quadrupole_angle
        self.MAX_angle_planarity = MAX_angle_planarity
        self.MAX_angle_Tshaped = MAX_angle_Tshaped

        
        #===== Create dictionaries of aromatic ring and corresponding plans =====
        self.dict_aromatic_ring = {"TYR": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
                                   "TRP": ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
                                   "PHE": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
                                   "HIS": ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
                                   "HID": ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
                                   "HIE": ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
                                   "HSD": ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
                                   "HSE": ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
                                  }
           
        self.dict_plane = {"TYR": ['CG',  'CE1', 'CE2'],
                           "TRP": ['CD2', 'CZ2', 'CZ3'],
                           "PHE": ['CG',  'CE1', 'CE2'],
                           "HIS": ['CG',  'CE1', 'NE2'],
                           "HID": ['CG',  'CE1', 'NE2'],
                           "HIE": ['CG',  'CE1', 'NE2'],
                           "HSE": ['CG',  'CE1', 'NE2'],
                           "HSD": ['CG',  'CE1', 'NE2'],
                          }
        
        
        
        #===== Identify residues properties and create their plan =====
        # get residue names
        self.res_index_A_name = self.top.residue(self.res_index_A).name
        self.res_index_B_name = self.top.residue(self.res_index_B).name
        
        # get atoms in the rings
        self.aromatic_A_atoms_list = self.dict_aromatic_ring[self.res_index_A_name]
        self.aromatic_B_atoms_list = self.dict_aromatic_ring[self.res_index_B_name]
        
        # Get position of atom makin the plane of the two residues
        self.list_aromatic_A_atom_position = []
        for self.atom in self.dict_plane[self.res_index_A_name]:
            self.atom_index = self.top.select(f"resid {self.res_index_A} and name {self.atom}")[0]
            self.atom_position = self.traj.xyz[0][self.atom_index]
            self.list_aromatic_A_atom_position.append(self.atom_position)

        self.list_aromatic_B_atom_position = []
        for self.atom in self.dict_plane[self.res_index_B_name]:
            self.atom_index = self.top.select(f"resid {self.res_index_B} and name {self.atom}")[0]
            self.atom_position = self.traj.xyz[0][self.atom_index]
            self.list_aromatic_B_atom_position.append(self.atom_position)    
            
        # Create the plane of each residue
        self.plane_aromatic_A = Plane.from_points(self.list_aromatic_A_atom_position[0], self.list_aromatic_A_atom_position[1], self.list_aromatic_A_atom_position[2])
        self.plane_aromatic_B = Plane.from_points(self.list_aromatic_B_atom_position[0], self.list_aromatic_B_atom_position[1], self.list_aromatic_B_atom_position[2])
        
        # Get normal vector of the planes
        self.vector_normal_plane_aromatic_A = self.plane_aromatic_A.normal
        self.vector_normal_plane_aromatic_B = self.plane_aromatic_B.normal
        
        
        #===== Distance COM-COM of aromatics =====
        # Compute COM of aromatic rings
        self.COM_aromatic_A = md.compute_center_of_mass(self.traj, select=f"resid {self.res_index_A} and name {' '.join(self.aromatic_A_atoms_list)}")[0]
        self.COM_aromatic_B = md.compute_center_of_mass(self.traj, select=f"resid {self.res_index_B} and name {' '.join(self.aromatic_B_atoms_list)}")[0]
        
        # create a vector between aromatic COMs and measure it's norm
        self.vector_COM_COM = Vector.from_points(self.COM_aromatic_A, self.COM_aromatic_B)
        self.distance_COM_COM = self.vector_COM_COM.norm() *10 # *10 to convert nm to angstrom
        
        
        
        #===== Distance between COM of resiaue A and the projected point of the COM of residue B into the plane of residue A =====
        # Get projected point of COM resid B on the plane of residue A
        self.projected_COM_aromatic_B = self.plane_aromatic_A.project_point(self.COM_aromatic_B)
        
        # create a vector between aromatic COM of residue A and the projected point
        self.vector_COM_projected = Vector.from_points(self.COM_aromatic_A, self.projected_COM_aromatic_B)
        self.distance_COM_projected = self.vector_COM_projected.norm() *10 # *10 to convert nm to angstrom

        
        
        #===== Angles (3) between norms and vector COM-COM =====
        
        #---- angle between the two plan normals -----
        # calculate the angle and convert it to degree
        self.angle_norm_norm = np.rad2deg( self.vector_normal_plane_aromatic_A.angle_between(self.vector_normal_plane_aromatic_B) )
        # Convert the angle to a range of values from 0 to 90 degrees
        self.angle_norm_norm = self.convert_angle(self.angle_norm_norm)
        
        #---- angle between the vector_COM_COM and normal of the plane of residue A -----
        # calculate the angle and convert it to degree
        self.angle_norm_A_COMCOM = np.rad2deg( self.vector_normal_plane_aromatic_A.angle_between(self.vector_COM_COM) )
        # Convert the angle to a range of values from 0 to 90 degrees
        self.angle_norm_A_COMCOM = self.convert_angle(self.angle_norm_A_COMCOM)      
        
        #---- angle between the vector_COM_COM and normal of the plane of residue B -----
        # calculate the angle and convert it to degree
        self.angle_norm_B_COMCOM = np.rad2deg( self.vector_normal_plane_aromatic_B.angle_between(self.vector_COM_COM) )
        # Convert the angle to a range of values from 0 to 90 degrees
        self.angle_norm_B_COMCOM = self.convert_angle(self.angle_norm_B_COMCOM)         
    
        #---- Convert angle to correspond to the plane and not the normal
        # The angle is before is between the normal vector of the plan and the vector COM --> charge
        # But the normal vector is perpendicular to the plan and geometric criteria are for the angle between the plan and the vector COM --> charge.
        #
        #    Pn                Pn is the normal vector to Plane.
        #    |   v             We calculate the angle between Pn and v.
        #    |  /              We need the angle between V and Plane. 
        #    | /
        # ___|/_____ Plane
        #
        # So angle is converted to correspond to be the one betwen the plan and the vector COM --> charge.
        self.angle_plane_plane = self.angle_norm_norm # angle between normals is the angle between the planes
        self.angle_plane_A_COMCOM = 90 - self.angle_norm_A_COMCOM
        self.angle_plane_B_COMCOM = 90 - self.angle_norm_B_COMCOM
    
    
        #===== Angles of Y or T shapes position =====
        #----- Generate all vector between the COM and the atoms in ring -----
        self.list_vector_COM_atom_ring_A = self.generate_list_vector_COM_atom_ring(self.COM_aromatic_A, self.res_index_A, self.aromatic_A_atoms_list)
        self.list_vector_COM_atom_ring_B = self.generate_list_vector_COM_atom_ring(self.COM_aromatic_B, self.res_index_B, self.aromatic_B_atoms_list)
    
        #----- Generate all angle between the normal of (other) plan and all vectors COM-atom_ring -----
        # Residue A
        self.list_angles_normal_vector_COM_atom_ring_A = []
        for self.vector in self.list_vector_COM_atom_ring_A:
            # calculate the angle and convert it to degree
            self.angle = np.rad2deg( self.vector_normal_plane_aromatic_B.angle_between(self.vector) )
            # Convert the angle to a range of values from 0 to 90 degrees
            self.angle = self.convert_angle(self.angle)
            # Convert angle to correspond to the plane and not the normal
            self.angle = 90 - self.angle
            # append the list
            self.list_angles_normal_vector_COM_atom_ring_A.append(self.angle)
         
        # Residue B
        self.list_angles_normal_vector_COM_atom_ring_B = []
        for self.vector in self.list_vector_COM_atom_ring_B:
            # calculate the angle and convert it to degree
            self.angle = np.rad2deg( self.vector_normal_plane_aromatic_A.angle_between(self.vector) )
            # Convert the angle to a range of values from 0 to 90 degrees
            self.angle = self.convert_angle(self.angle)
            # Convert angle to correspond to the plane and not the normal
            self.angle = 90 - self.angle
            # append the list
            self.list_angles_normal_vector_COM_atom_ring_B.append(self.angle)
        
        # Get the smallest angle from 'self.list_angles_normal_vector_COM_atom_ring_A' and 'self.list_angles_normal_vector_COM_atom_ring_B'
        self.angle_shaped = min(self.list_angles_normal_vector_COM_atom_ring_A + self.list_angles_normal_vector_COM_atom_ring_B)
    
    
    
    #===== Functions =====
    def convert_angle(self, angle):
        """
        Take angle in degree, in the range 0-360, and convert it to an equavalent range 0-90.
        """
        if angle <= 90:
            return angle # the angle is already in the range 0 to 90 degrees
        elif angle <= 180:
            return 180 - angle
        elif angle <= 270:
            return angle - 180
        else:
            return 360 - angle
    
    #-----    
    def generate_list_vector_COM_atom_ring(self, COM_position, residue_ID, aromatic_atoms_list):
        """
        Create a list of vector COM-atom for all atom in the aromatic ring
        """
        # create a list to store all vectors
        list_vector_ring = []
        
        # generate all possible COM-atom_ring vectors
        for atom in aromatic_atoms_list:
            # Get position of the atom position
            atom_index = self.top.select(f"resid {residue_ID} and name {atom}")[0]
            atom_position = self.traj.xyz[0][atom_index]
            
            # create the vector COM-atom
            vector = Vector.from_points(COM_position, atom_position)
            
            # Add the vector to the list
            list_vector_ring.append(vector)
            
        # return the list containing all vectors
        return list_vector_ring
    
    
    
    #===== Return properties =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if two aromatic residues interact with each other.
        
        RETURN  
            True, 'subtype'  The interaction exist            
            False, False     The interaction don't exist.
        """
        # check if the distance between the COM-COM is greeter than the 'MAX_distance_COM'
        if self.distance_COM_COM > self.MAX_distance_COM:
            return False, False
        
        # if the distance is not greter than 'MAX_distance_COM', the identify the subtype
        else:
            
            # identify coplanar parameters
            if self.angle_plane_plane <= self.MAX_angle_planarity and self.angle_plane_A_COMCOM <= self.MAX_quadrupole_angle and self.angle_plane_B_COMCOM <= self.MAX_quadrupole_angle:
                return True, 'coplanar'
            
             # identify parallel parameters
            elif self.angle_plane_plane <= self.MAX_angle_planarity and self.MIN_pi_angle <= self.angle_plane_A_COMCOM and self.MIN_pi_angle <= self.angle_plane_B_COMCOM \
            and self.distance_COM_projected < self.MIN_distance_offset:
                return True, 'parallel'
            
             # identify offset parameters    
            elif self.angle_plane_plane <= self.MAX_angle_planarity and self.MIN_pi_angle <= self.angle_plane_A_COMCOM and self.MIN_pi_angle <= self.angle_plane_B_COMCOM \
            and self.MIN_distance_offset <= self.distance_COM_projected <= self.MAX_distance_offset:
                return True, 'offset'
            
            # identify T/Y-shaped parameters
            elif self.MIN_pi_angle <= self.angle_plane_plane \
            and ((self.MIN_pi_angle <= self.angle_plane_A_COMCOM and self.angle_plane_B_COMCOM <= self.MAX_quadrupole_angle) or (self.MIN_pi_angle <= self.angle_plane_B_COMCOM and self.angle_plane_A_COMCOM <= self.MAX_quadrupole_angle)):
                
                # the type is T-shaped or Y-shaped. Now distinguish between them:
                
                # if at least one of the angle betwwen vector COM-atom_ring and the opposit ring plane correspond to the pi area
                if self.angle_shaped <= self.MAX_angle_Tshaped:
                    return True, 'T-shaped'
                else:
                    return True, 'Y-shaped'
                
            # identify intermediate parameters
            else:
                return True, 'intermediate'
        



    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         [atoms_interaction]
        """
        atoms_interaction = []

        aromatic_ring_A = list(self.top.select(f"resid {self.res_index_A} and name {' '.join(self.aromatic_A_atoms_list)}"))
        aromatic_ring_B = list(self.top.select(f"resid {self.res_index_B} and name {' '.join(self.aromatic_B_atoms_list)}"))

        atoms_interaction.extend(aromatic_ring_A)
        atoms_interaction.extend(aromatic_ring_B)

        return [atoms_interaction]



    @property
    def get_angle(self):
        """
        DESCRIPTION    Angle between the atomatic plane, and the vector COM-COM with each of the plane
        RETURN         angle_plane_plane, angle_plane_A_COMCOM, angle_plane_B_COMCOM, angle_shaped
        UNIT           degree
        """
        return self.angle_plane_plane, self.angle_plane_A_COMCOM, self.angle_plane_B_COMCOM, self.angle_shaped

    
    @property
    def get_distance(self):
        """
        DESCRIPTION    COM-COM distance between the aromatic rings, and their COM-projectedCOM distance
        RETURN         distance_COM_COM, distance_COM_projected
        UNIT           Angstrom
        """
        return self.distance_COM_COM, self.distance_COM_projected

    



#=====================================================
#===== Class for Arg-Arg & Arg-Aromatic
#=====================================================    
class arg_involved:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_distance=6.0, MIN_perpendicular_angle=60.0, MAX_parallel_angle=35.0):
        """
        INTERACTION TYPE    Arginine involved stacking: Arg-Aromatic, Arg-Arg
        SUBTYPE(S)          perpendicular, parallel, intermediate

        DESCRIPTION
            Interaction with the charged amino acid and an aromatic rings.
            Please note that protonated histidine are not taken in acount are not taken into account because they can make charge-aromatic interactions.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS
            frame    Frame ID on which to perform the analysis.
                     Default value: 0
            
            MAX_distance    Distance between the CZ atom of the ARG and COM of and aromatic ring or CZ of another ARG.
                            Default value: 6.0 Å

            MIN_perpendicular_angle    Minimum angle defining if the planes are perpendicular or not (the maximum is 90˚).
                                       Default value: 60.0˚

            MAX_parallel_angle    Maximum angle defining if the planes are parallel or not (the maximum is 0˚).
                                  Default value: 35.0˚
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.MAX_distance = MAX_distance
        self.MIN_perpendicular_angle = MIN_perpendicular_angle
        self.MAX_parallel_angle = MAX_parallel_angle
        self.res_index_A = res_index_A
        self.res_index_B = res_index_B

        
        #===== Create dictionaries of aromatic ring and corresponding plans =====
        self.dict_aromatic_ring = {"TYR": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "TRP": 'CD2 CE2 CE3 CZ2 CZ3 CH2',
                                   "PHE": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "HIS": 'CG ND1 CD2 CE1 NE2',
                                   "HID": 'CG ND1 CD2 CE1 NE2',
                                   "HIE": 'CG ND1 CD2 CE1 NE2',
                                   "HSD": 'CG ND1 CD2 CE1 NE2',
                                   "HSE": 'CG ND1 CD2 CE1 NE2',
                                  }
           
        self.dict_plane = {"TYR": ['CG',  'CE1', 'CE2'],
                           "TRP": ['CD2', 'CZ2', 'CZ3'],
                           "PHE": ['CG',  'CE1', 'CE2'],
                           "HIS": ['CG',  'CE1', 'NE2'],
                           "HID": ['CG',  'CE1', 'NE2'],
                           "HIE": ['CG',  'CE1', 'NE2'],
                           "HSE": ['CG',  'CE1', 'NE2'],
                           "HSD": ['CG',  'CE1', 'NE2'],
                           "ARG": ['NE',  'NH1', 'NH2'],
                          }
        
        
        
        #===== Identify residues properties and create the type of interaction =====       
        if self.top.residue(res_index_A).name == 'ARG' and self.top.residue(res_index_B).name in self.dict_aromatic_ring:
            self.res_arg_index = res_index_A
            self.res_aromatic_index = res_index_B
            self.res_aromatic_name = self.top.residue(res_index_B).name
            self.type = "Arg-Aromatic"
            
        elif self.top.residue(res_index_B).name == 'ARG' and self.top.residue(res_index_A).name in self.dict_aromatic_ring:
            self.res_arg_index = res_index_B
            self.res_aromatic_index = res_index_A
            self.res_aromatic_name = self.top.residue(res_index_A).name
            self.type = "Arg-Aromatic"
            
        elif self.top.residue(res_index_A).name == 'ARG' and self.top.residue(res_index_B).name == 'ARG':
            self.res_index_A = res_index_A
            self.res_index_B = res_index_B
            self.type = "Arg-Arg"
            
        else:
            raise ValueError("Residue A and B must be both ARG, or one must be ARG and the other is aromatic.")
        
        
        
        #===== Identify interaction for "Arg-Aromatic"  =====
        if self.type == "Arg-Aromatic":
            
            #----- Get angle between ARG plane and aromatic plane -----
            # get atoms in the aromatic ring
            #self.aromatic_atoms_list = self.dict_aromatic_ring[self.res_aromatic_name]
            
            # Get position of atom making the aromatic plane
            self.list_aromatic_atom_position = []
            for self.atom in self.dict_plane[self.res_aromatic_name]:
                self.atom_index = self.top.select(f"resid {self.res_aromatic_index} and name {self.atom}")[0]
                self.atom_position = self.traj.xyz[0][self.atom_index]
                self.list_aromatic_atom_position.append(self.atom_position)
                
            # Get position of atom making the ARG plane
            self.list_arg_atom_position = []
            for self.atom in self.dict_plane['ARG']:
                self.atom_index = self.top.select(f"resid {self.res_arg_index} and name {self.atom}")[0]
                self.atom_position = self.traj.xyz[0][self.atom_index]
                self.list_arg_atom_position.append(self.atom_position)
            
            # Create the 2 planes
            self.plane_aromatic = Plane.from_points(self.list_aromatic_atom_position[0], self.list_aromatic_atom_position[1], self.list_aromatic_atom_position[2])
            self.plane_arg = Plane.from_points(self.list_arg_atom_position[0], self.list_arg_atom_position[1], self.list_arg_atom_position[2])
            
            # Get normal vector of the planes
            self.vector_normal_plane_aromatic = self.plane_aromatic.normal
            self.vector_normal_plane_arg = self.plane_arg.normal
            
            # Calculate the angle between the two planes
            self.angle = np.rad2deg( self.vector_normal_plane_aromatic.angle_between(self.vector_normal_plane_arg) )
            # Convert the angle to a range of values from 0 to 90 degrees
            self.angle = self.convert_angle(self.angle)
            
            #----- Measure distance between CZ of ARG and COM of aromatic ring -----
            # Get 
            self.COM_aromatic = md.compute_center_of_mass(self.traj, select=f"resid {self.res_aromatic_index} and name {self.dict_aromatic_ring[self.res_aromatic_name]}")[0]
            self.CZ_atom_arg_index = self.top.select(f"resid {self.res_arg_index} and name CZ")[0]
            self.CZ_atom_arg_position = self.traj.xyz[0][self.CZ_atom_arg_index]
            
            # Get the distance between the aromatic ring COM and the CZ atom of ARG
            self.vector_COM_CZ = Vector.from_points(self.COM_aromatic, self.CZ_atom_arg_position)
            self.distance = self.vector_COM_CZ.norm() *10 # *10 to convert nm to angstrom
        
            
        
        #===== Identify interaction for "Arg-Arg"  =====
        elif self.type == "Arg-Arg":
            #----- Calculate angle between the 2 planes -----
            # Get position of atom making the ARG plane A
            self.list_arg_A_atom_position = []
            for self.atom in self.dict_plane['ARG']:
                self.atom_index = self.top.select(f"resid {self.res_index_A} and name {self.atom}")[0]
                self.atom_position = self.traj.xyz[0][self.atom_index]
                self.list_arg_A_atom_position.append(self.atom_position)
                
            # Get position of atom making the ARG plane B
            self.list_arg_B_atom_position = []
            for self.atom in self.dict_plane['ARG']:
                self.atom_index = self.top.select(f"resid {self.res_index_B} and name {self.atom}")[0]
                self.atom_position = self.traj.xyz[0][self.atom_index]
                self.list_arg_B_atom_position.append(self.atom_position)
 
            # Create the 2 planes
            self.plane_A = Plane.from_points(self.list_arg_A_atom_position[0], self.list_arg_A_atom_position[1], self.list_arg_A_atom_position[2])
            self.plane_B = Plane.from_points(self.list_arg_B_atom_position[0], self.list_arg_B_atom_position[1], self.list_arg_B_atom_position[2])
            
            # Get normal vector of the planes
            self.vector_normal_plane_A = self.plane_A.normal
            self.vector_normal_plane_B = self.plane_B.normal
            
            # Calculate the angle between the two planes
            self.angle = np.rad2deg( self.vector_normal_plane_A.angle_between(self.vector_normal_plane_B) )
            # Convert the angle to a range of values from 0 to 90 degrees
            self.angle = self.convert_angle(self.angle)
            
            #----- Calculate the distance between the two ARG -----
            # Get 
            self.CZ_A_index = self.top.select(f"resid {self.res_index_A} and name CZ")[0]
            self.CZ_A_position = self.traj.xyz[0][self.CZ_A_index]
            self.CZ_B_index = self.top.select(f"resid {self.res_index_B} and name CZ")[0]
            self.CZ_B_position = self.traj.xyz[0][self.CZ_B_index]
            
            # Get the distance between the aromatic ring COM and the CZ atom of ARG
            self.vector_CZ = Vector.from_points(self.CZ_A_position, self.CZ_B_position)
            self.distance = self.vector_CZ.norm() *10 # *10 to convert nm to angstrom            
    
    
    #===== Functions =====
    def convert_angle(self, angle):
        """
        Take angle in degree, in the range 0-360, and convert it to an equavalent range 0-90.
        """
        if angle <= 90:
            return angle # the angle is already in the range 0 to 90 degrees
        elif angle <= 180:
            return 180 - angle
        elif angle <= 270:
            return angle - 180
        else:
            return 360 - angle
    
    
    
    #===== Return properties =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if a charged residue interact between Arg-Arg or Arg-Aromatic
        
        RETURN  
            True, 'subtype'  The interaction exist            
            False, False     The interaction don't exist.
        """
        if self.distance > self.MAX_distance:
            return False, False, False
        
        else:
            if self.angle <= self.MAX_parallel_angle:
                return True, self.type, "parallel"
            
            elif self.MIN_perpendicular_angle <= self.angle:
                return True, self.type, "perpendicular"
            
            else:
                return True, self.type, "intermediate"

        
    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         [atoms_interaction] or [[CZ_A_index, CZ_B_index]]
        """
        if self.type == "Arg-Aromatic":
            atoms_interaction = list(self.top.select(f"resid {self.res_aromatic_index} and name {self.dict_aromatic_ring[self.res_aromatic_name]}"))
            atoms_interaction.append(self.CZ_atom_arg_index)
            return [atoms_interaction]

        elif self.type == "Arg-Arg":
            return [[self.CZ_A_index, self.CZ_B_index]]




    @property
    def get_angle(self):
        """
        DESCRIPTION    Angle between the 2 Arg plane or the angle between the Arg plane and the aromatic plane
        RETURN         angle
        UNIT           degree
        """
        return self.angle

    
    @property
    def get_distance(self):
        """
        DESCRIPTION    Distance between the CZ of each Arg, of distance between the CZ of Arg and the COM of the aromatic plane
        RETURN         distance
        UNIT           Angstrom
        """
        return self.distance
    
    




#=====================================================
#===== Class for Pi-Hbond
#=====================================================
class pi_hbond:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MIN_angle=120.0, MAX_distance_X_COM=5.5, MAX_distance_H_COM=3.0, MAX_distance_Hp_COM=1.2):
        """
        INTERACTION TYPE    Stacking of two aromatic residues.
        SUBTYPE(S)          

        DESCRIPTION
            Stacking of two aromatic rings.
            Please note that protonated histidine are not taken in acount are not taken into account because they are involved in charge-aromatic interactions.

            * COM is the center  of mass of the aromatic ring 
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS
            frame    Frame ID on which to perform the analysis.
                     Default value: 0

            MIN_angle    Minimum angle of X-H...COM
                         Default value: 120.0
            
            MAX_distance_X_COM    Maximum distance between H donnor (X) and COM
                                  Default value: 5.5

            MAX_distance_H_COM    Maximum distance between H and COM
                                  Default value: 3.0

            MAX_distance_Hp_COM   Maximum distance between H projected on aromatic plane and COM
                                  Default value: 1.2
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.MIN_angle = MIN_angle
        self.MAX_distance_X_COM = MAX_distance_X_COM
        self.MAX_distance_H_COM = MAX_distance_H_COM
        self.MAX_distance_Hp_COM = MAX_distance_Hp_COM 

        
        #===== Create dictionaries of aromatic ring and corresponding plans =====
        #----- C carbon in aromatic ring -----
        self.dict_aromatic_ring = {"TYR": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "TRP": 'CD2 CE2 CE3 CZ2 CZ3 CH2',
                                   "PHE": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "HIS": 'CG ND1 CD2 CE1 NE2',
                                   "HID": 'CG ND1 CD2 CE1 NE2',
                                   "HIE": 'CG ND1 CD2 CE1 NE2',
                                   "HSD": 'CG ND1 CD2 CE1 NE2',
                                   "HSE": 'CG ND1 CD2 CE1 NE2',
                                  }
        
                   
        self.dict_plane = {"TYR": ['CG',  'CE1', 'CE2'],
                           "TRP": ['CD2', 'CZ2', 'CZ3'],
                           "PHE": ['CG',  'CE1', 'CE2'],
                           "HIS": ['CG',  'CE1', 'NE2'],
                           "HID": ['CG',  'CE1', 'NE2'],
                           "HIE": ['CG',  'CE1', 'NE2'],
                           "HSE": ['CG',  'CE1', 'NE2'],
                           "HSD": ['CG',  'CE1', 'NE2'],
                          }
        
        #----- Ignored H: H in aromatic ring (due to aromatic-aromatic stacking how can involve H) & H in charged residue (due to the charge-aromatic interaction how can involve H) -----
        self.dict_ignored_H = {"TYR": "HE1 HE2 HD1 HD2 HH",
                               "TRP": "ND1 HE1 HZ2 HH2 HZ3 HE3",
                               "PHE": "HE1 HE2 HD1 HD2 HZ",
                               "HIS": "HE1 HE2 HD1 HD2",
                               "HID": "HE1 HE2 HD1 HD2",
                               "HIE": "HE1 HE2 HD1 HD2",
                               "HSD": "HE1 HE2 HD1 HD2",
                               "HSE": "HE1 HE2 HD1 HD2",
                               "HIP": "HE1 HE2 HD1 HD2", # Protonated HIS
                               "HSP": "HE1 HE2 HD1 HD2", # Protonated HIS in CHARMM
                               "ARG": "HE HH11 HH12 HH21 HH22",
                               "LYS": "HZ1 HZ2 HZ3",
                               "ASP": "HD2",
                               "GLU": "HE2",
                              }
 

        #===== initialise residue name and index =====
        if self.top.residue(res_index_A).name not in self.dict_aromatic_ring and self.top.residue(res_index_B).name in self.dict_aromatic_ring:
            self.res_R_index = res_index_A
            self.res_R_name = self.top.residue(res_index_A).name
            self.res_aromatic_index = res_index_B
            self.res_aromatic_name = self.top.residue(res_index_B).name
        
        elif self.top.residue(res_index_B).name not in self.dict_aromatic_ring and self.top.residue(res_index_A).name in self.dict_aromatic_ring:
            self.res_R_index = res_index_B
            self.res_R_name = self.top.residue(res_index_B).name
            self.res_aromatic_index = res_index_A
            self.res_aromatic_name = self.top.residue(res_index_A).name 

        
        #===== COM and plane of aromatic ring =====
        # Calculate the COM of the aromatic ring
        self.COM_aromatic = md.compute_center_of_mass(self.traj, select=f"resid {self.res_aromatic_index} and name {self.dict_aromatic_ring[self.res_aromatic_name]}")[0]
        
        # Get position of atom making the aromatic plane
        self.list_aromatic_atom_position = []
        for self.atom in self.dict_plane[self.res_aromatic_name]:
            self.atom_index = self.top.select(f"resid {self.res_aromatic_index} and name {self.atom}")[0]
            self.atom_position = self.traj.xyz[0][self.atom_index]
            self.list_aromatic_atom_position.append(self.atom_position)
        
        # Create the plane
        self.plane_aromatic = Plane.from_points(self.list_aromatic_atom_position[0], self.list_aromatic_atom_position[1], self.list_aromatic_atom_position[2])
        
        
        
        #===== Identify all Hbonds =====
        self.checking = False
        
        #----- Get all H indices of X residue -----
        if self.res_R_name not in self.dict_ignored_H:
            self.hydrogen_indices = self.top.select(f"resid {self.res_R_index} and element H")
        else:
            self.hydrogen_indices = self.top.select(f"resid {self.res_R_index} and element H and not name {self.dict_ignored_H[self.res_R_name]}")
        
        #----- create list to store results -----
        self.list_XH_indices = []
        self.list_distances = []
        self.list_angles = []
        
        #----- Identify posibles Hbonds using distances and angles -----
        for self.H_index in self.hydrogen_indices:
            #+++ 1. get position of the H atom  and it's atom in the topology
            self.H_position = self.traj.xyz[0][self.H_index]
            self.atom_H = self.traj.topology.atom(self.H_index)
            
            
            #+++ 2. get distances COM-H and COM-Hp (Hp: H projected on the aromatic plane)
            # create a vector between the COM and the H, and get it's norm
            self.vector_H_COM = Vector.from_points(self.H_position, self.COM_aromatic)
            self.distance_H_COM = self.vector_H_COM.norm() *10 # *10 to convert nm to angstrom
                
            # project H on the aromatic plane
            self.H_projected_position = self.plane_aromatic.project_point(self.H_position)
            
            # create a vector between the COM and the H projected on the plane, and get it's norm
            self.vector_COM_Hprojected = Vector.from_points(self.COM_aromatic, self.H_projected_position)
            self.distance_COM_Hprojected = self.vector_COM_Hprojected.norm() *10 # *10 to convert nm to angstrom  
            
            
            #+++ 3. identify X atom bonded to H and the distance X-COM
            # get the covalent bond involving H
            for self.bond in self.top.bonds:
                if self.atom_H in self.bond:
                    self.bond_XH = self.bond
            
            if self.bond_XH[0] != self.atom_H:
                self.atom_X = self.bond_XH[0]
            else:
                self.atom_X = self.bond_XH[1]
            
            # get the position of the X atom
            self.X_index = self.atom_X.index
            self.X_position = self.traj.xyz[0][self.X_index]
            
            # Get the distance X-COM
            self.vector_X_COM = Vector.from_points(self.X_position, self.COM_aromatic)
            self.distance_X_COM = self.vector_X_COM.norm() *10 # *10 to convert nm to angstrom
            
            
            #+++ 4. calculate the X-H-COM angle
            self.vector_H_X = Vector.from_points(self.H_position, self.X_position)
            self.angle = np.rad2deg( self.vector_H_X.angle_between(self.vector_H_COM) )
            
            
            #+++ 5. select the bond using COM-H, COM-Hp, COM-X distances and the angle X-H-COM           
            self.list_XH_indices.append([self.H_index, self.X_index])
            self.list_distances.append([self.distance_H_COM, self.distance_COM_Hprojected, self.distance_X_COM])
            self.list_angles.append([self.angle])
            
            if self.distance_H_COM <= self.MAX_distance_H_COM and self.distance_COM_Hprojected <= self.MAX_distance_Hp_COM \
            and self.distance_X_COM <= self.MAX_distance_X_COM and self.MIN_angle <= self.angle:
                self.checking = True
                
    
    
    #===== Return results =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the amino group of GLN or ASN interact with an aromatic ring.
        
        RETURN
            True     The interaction exist.
            False    The interaction don't exist.
        """
        return self.checking

    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return a list of atoms index involved in H-bond.
                       X atom is the hydrogen donor.
        RETURN         [[H_index, X_index],...]
        """
        return self.list_XH_indices        
        
    @property
    def get_angle(self):
        """
        DESCRIPTION    Angle betwee the N-COM and the normal of the aromatic ring
        RETURN         [[angle],...]
        UNIT           degree
        """
        return self.list_angles

    @property
    def get_distance(self):
        """
        DESCRIPTION    Distance between COM of aromatic ring and N of the amino group
        RETURN         [[distance_H_COM, distance_COM_Hprojected, distance_X_COM],...]
        UNIT           Angstrom
        """
        return self.list_distances           
    




#=====================================================
#===== Class for n-->pi*
#=====================================================
class n_pi:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, ref_distance=3.0, distance_tolerance=0.25, ref_angle=110.0 , angular_tolerance=10.0):
        """
        INTERACTION TYPE    n --> Pi*
        SUBTYPE(S)          regular or reciprocal

        DESCRIPTION
            Interaction (C=O...C=O) between tow carbonyl groups
            (C=O) of the backbone.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS
        Set an absolute tolerance parameter N. See documentation concerning. 'numpy.isclose'.

            ref_distance          Distance between O..C atoms.
                                  Default value: 3.0 Å

            distance_tolerance    The range of distance is ref_distance +/- N.
                                  Default value of N: 0.25

            ref_angle             Angle between O..C=O atoms.
                                  Default value: 110.0

            angular_tolerance     The range of angle is ref_angle +/- N.
                                  Default value of N: 10.0

            frame                 Frame ID on which to perform the analysis.
                                  Default value: 0
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        self.ref_distance = ref_distance
        self.distance_tolerance = distance_tolerance
        self.angular_tolerance = angular_tolerance
        self.ref_angle = ref_angle

        #===== Get atoms index for residue A =====
        self.atom_O_res_A = self.top.select(f"resid {self.res_A} and name O")[0]
        self.atom_C_res_A = self.top.select(f"resid {self.res_A} and name C")[0]
        
        #===== Get atoms index for residue B =====
        self.atom_O_res_B = self.top.select(f"resid {self.res_B} and name O")[0]
        self.atom_C_res_B = self.top.select(f"resid {self.res_B} and name C")[0]
        
        #===== Compute O...C distances =====
        # *10 is used to convert nm to angstrom
        self.distance_O_res_A_C_res_B = md.compute_distances(self.traj, [[self.atom_O_res_A, self.atom_C_res_B]])[0][0] *10
        self.distance_O_res_B_C_res_A = md.compute_distances(self.traj, [[self.atom_O_res_B, self.atom_C_res_A]])[0][0] *10
        
        #===== Compute O...C=O angles =====
        self.angle_O_res_A_CO_res_B = np.rad2deg( md.compute_angles(self.traj, [[self.atom_O_res_A, self.atom_C_res_B, self.atom_O_res_B]])[0][0] )
        self.angle_O_res_B_CO_res_A = np.rad2deg( md.compute_angles(self.traj, [[self.atom_O_res_B, self.atom_O_res_A, self.atom_C_res_A]])[0][0] )
    
    
    
    #===== Return results =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the n --> Pi* interaction exist between tow residues.
        
        RETURN
            If the interaction exist between the 2 residues:
                True, reregular
                True, reciprocal
            
            If the interaction don't exist between the 2 residues:
                False, False
        """
        
        # If the two distances are close to self.ref_distance and the two angles are close to self.ref_angle
        # The interaction is reciprocal
        if (np.isclose(self.distance_O_res_A_C_res_B, self.ref_distance, atol=self.distance_tolerance) and np.isclose(self.angle_O_res_A_CO_res_B, self.ref_angle, atol=self.angular_tolerance)) \
        and (np.isclose(self.distance_O_res_B_C_res_A, self.ref_distance, atol=self.distance_tolerance) and np.isclose(self.angle_O_res_B_CO_res_A, self.ref_angle, atol=self.angular_tolerance)):
            return True, 'reciprocal'
        
        # # If only one residue have a distance close to self.ref_distance and an angle close to 102.0
        # The interaction is regular
        elif (np.isclose(self.distance_O_res_A_C_res_B, self.ref_distance, atol=self.distance_tolerance) and np.isclose(self.angle_O_res_A_CO_res_B, self.ref_angle, atol=self.angular_tolerance)) \
        or (np.isclose(self.distance_O_res_B_C_res_A, self.ref_distance, atol=self.distance_tolerance) and np.isclose(self.angle_O_res_B_CO_res_A, self.ref_angle, atol=self.angular_tolerance)):
            return True, 'regular'
        
        # If the two distances are different to 3.0 and the two angles are different to 102.0
        # The interaction don't exist.
        else:
            return False, False
        


    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         [[atom_O_res_A,atom_C_res_B],[atom_O_res_B,atom_C_res_A]]
        """
        return [[self.atom_O_res_A,self.atom_C_res_B],[self.atom_O_res_B,self.atom_C_res_A]]
    
    
    
    #===== Return geometrical properties =====
    @property
    def get_distance(self):
        """
        DESCRIPTION    Return the two O...C distances.
        RETURN         distance_O_res_A_C_res_B, distance_O_res_B_C_res_A
        """
        return self.distance_O_res_A_C_res_B, self.distance_O_res_B_C_res_A
    
    
    @property
    def get_angle(self):
        """
        DESCRIPTION    Return the two O...C=O angles.
        RETURN         angle_O_res_A_CO_res_B, angle_O_res_B_CO_res_A
        """
        return self.angle_O_res_A_CO_res_B, self.angle_O_res_B_CO_res_A
    




#=====================================================
#===== Class for Chalcogen bond & S/Se mediated H-bond
#=====================================================
class sse_hydrogen_chalcogen_bond:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_distance_SX=3.6, MIN_angle_theta_chalcogen = 50.0, MAX_angle_dihedral_chalcogen=50.0, \
                  MIN_angle_phi_chalcogen=30.0, MAX_angle_phi_chalcogen=60.0, MIN_angle_csx_chalcogen=115.0, MAX_angle_csx_chalcogen=155.0, MAX_angle_csx_hbond=145.0):
        """
        INTERACTION TYPE    Chalcogen bond or S/Se mediated H-bond
        SUBTYPE(S)          no

        DESCRIPTION
            Chalogen bond: sigma-hole interaction between S/Se and N/O.
            
            S/Se mediated H-bond: Hydrogen bond for which S/Se is acceptor or donnor.
       
        REMARKS
            1) Analysis don't take in account the presence of hydrogen.
            2) 'IndexError: list index out of range' Always happend when bonded atoms are missing.
               Ensure to have CONNECT section in your PDB, with the correct bonds.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS
            frame                 Frame ID on which to perform the analysis.
                                  Default value: 0
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.res_A = res_index_A
        self.res_B = res_index_B
        self.MAX_distance_SX = MAX_distance_SX
        self.MAX_angle_dihedral_chalcogen = MAX_angle_dihedral_chalcogen
        self.MIN_angle_phi_chalcogen = MIN_angle_phi_chalcogen
        self.MAX_angle_phi_chalcogen = MAX_angle_phi_chalcogen
        self.MIN_angle_theta_chalcogen = MIN_angle_theta_chalcogen
        self.MIN_angle_csx_chalcogen = MIN_angle_csx_chalcogen
        self.MAX_angle_csx_chalcogen = MAX_angle_csx_chalcogen
        self.MAX_angle_csx_hbond = MAX_angle_csx_hbond
        
        #===== Create list to store results =====
        self.list_chalcogen = []
        self.list_chalcogen_distance = []
        self.list_chalcogen_angle = []
        self.list_hbond = []
        self.list_hbond_distance = []
        self.list_hbond_angle = []
        
        
        #===== Get all S/Se atoms =====
        self.list_index_atom_SSE = self.top.select(f"resid {self.res_A} {self.res_B} and element S Se")
        
        
        #===== For all S/Se atoms get geometric parameters =====
        #----- Get Points P1 and P2 to create the plane ----- 
        # 
        for self.index_atom_SSE in self.list_index_atom_SSE:
            
            # Get the atom S/SE in topology
            self.atom_SSE = self.top.atom(self.index_atom_SSE) 
            
            # create a list to store atom covalent to S/Se
            self.list_atom_bonded_to_SSE = []
            
            # Search bonded atoms to S/Se of the residue
            for self.bond in self.traj.topology.bonds:
                # Check if bonded to S/Se atom an append the list_atom_bonded_to_SSE
                if self.atom_SSE in self.bond:
                    if self.atom_SSE != self.bond[0]: #if bond[0] is not the S/Se atom, append the list
                        self.list_atom_bonded_to_SSE.append(self.bond[0])
                        
                    elif self.atom_SSE != self.bond[1]: #else if bond[1] is not the S/Se atom, append the list
                        self.list_atom_bonded_to_SSE.append(self.bond[1])                
            
            # remove duplicate in the list
            self.list_atom_bonded_to_SSE = list(set(self.list_atom_bonded_to_SSE))
        
        
            #----- Get atom coordinates -----
            self.coordinates_atom_P1 = self.traj.xyz[0, self.list_atom_bonded_to_SSE[0].index]
            self.coordinates_atom_SSE = self.traj.xyz[0, self.index_atom_SSE]
            self.coordinates_atom_P2 = self.traj.xyz[0, self.list_atom_bonded_to_SSE[1].index]
            
            #----- Get the centorid and make the plane of P1-SSE-P2 -----
            self.centroid = Points([self.coordinates_atom_P1, self.coordinates_atom_SSE, self.coordinates_atom_P2]).centroid()
            self.plane = Plane.from_points(self.coordinates_atom_SSE, self.coordinates_atom_P1, self.coordinates_atom_P2)
            
            #----- Create vectors -----
            self.vec_centroid_SSE = Vector.from_points(self.centroid, self.coordinates_atom_SSE)
            self.vec_SSE_centroid = Vector.from_points(self.coordinates_atom_SSE, self.centroid) # vector with opposit direction of vec_centroid_SSE
            
            #----- Get neighbors atom N,O (not H) -----
            # Get the residue (residue_opposit) that don't contain the atom_SSE
            if self.atom_SSE.residue.index != self.res_A:
                self.residue_opposit = self.top.residue(self.res_A)
            
            elif self.atom_SSE.residue.index != self.res_B:
                self.residue_opposit = self.top.residue(self.res_B)
            
            # For all N and O in residue_opposit, measure the distance with the atom S/Se
            self.list_atoms_close_to_SSE = []
            self.list_distances_atoms_close_to_SSE = []
            
            for self.atom_X in self.residue_opposit.atoms:
                if self.atom_X.element.symbol in ['N','O']:
                    self.distance = md.compute_distances(self.traj, [[self.index_atom_SSE, self.atom_X.index]])[0][0] *10
                    
                    # Keep only atoms with distance <= MAX_distance_SX to S/Se atom
                    if self.distance <= self.MAX_distance_SX:
                        
                        # Append results list with distance and atom index
                        self.list_atoms_close_to_SSE.append(self.atom_X.index)
                        self.list_distances_atoms_close_to_SSE.append(self.distance)
                        
                        # Get coordiante of atom_X and it's projection on the plane P1-SSE-P2
                        self.coordinates_atom_X = self.traj.xyz[0, self.atom_X.index]
                        self.coordinates_projection_atom_X = self.plane.project_point(self.coordinates_atom_X)
                        
                        # Creates vectors
                        self.vec_SSE_atom_X = Vector.from_points(self.coordinates_atom_SSE, self.coordinates_atom_X)
                        self.vec_SSE_projection_atom_X = Vector.from_points(self.coordinates_atom_SSE, self.coordinates_projection_atom_X)
                        
                        # Compute angles, and convert then into degreee, then in range 0-180 with the command self.convert_angle_180()
                        #..... angle between the normal of the plan and the vector SSE --> atom_X
                        self.angle_theta = self.convert_angle_180(np.rad2deg( self.plane.normal.angle_between(self.vec_SSE_atom_X) ))
                        
                        #..... angle between the vector centroid-->SSE and vector SSE-->projected_atom_X
                        self.angle_phi = self.convert_angle_180(np.rad2deg( self.vec_centroid_SSE.angle_between(self.vec_SSE_projection_atom_X) ))
                        
                        #..... angle between vectors centroid-->SSE and vector SSE-->atom_X
                        self.angle_CSX = self.convert_angle_180(np.rad2deg( self.vec_SSE_centroid.angle_between(self.vec_SSE_atom_X) ))
                       
                        #..... identify carbon atom (P) bonded to S/Se and make the plane P-centroid-SSE
                        #      and create vector SSE-->P
                        if self.list_atom_bonded_to_SSE[0].name in ['CB', 'CG']:
                            self.plane_P = Plane.from_points(self.coordinates_atom_P1, self.centroid, self.coordinates_atom_SSE)
                        
                        elif self.list_atom_bonded_to_SSE[1].name in ['CB', 'CG']:
                            self.plane_P = Plane.from_points(self.coordinates_atom_P2, self.centroid, self.coordinates_atom_SSE)
                    
                        #..... dihedral angle P--centroid--SSE--atom_X
                        # create plane: centroid-SSE-atom_X
                        self.plane_CSX = Plane.from_points(self.centroid, self.coordinates_atom_SSE, self.coordinates_atom_X)
                        
                        # compute dihedral angles = 180 - angle between the normals, then ensure they are in range 0-180
                        self.angle_dihedral = self.convert_angle_180(180 - np.rad2deg( self.plane_P.normal.angle_between( self.plane_CSX.normal )))
                        
                        # Check chalcogen bond
                        if self.angle_dihedral <= self.MAX_angle_dihedral_chalcogen and self.MIN_angle_theta_chalcogen < self.angle_theta \
                        and (self.MIN_angle_phi_chalcogen <= self.angle_phi <= self.MAX_angle_phi_chalcogen) \
                        and (self.MIN_angle_csx_chalcogen <= self.angle_CSX <= self.MAX_angle_csx_chalcogen):
                            
                            # append list of atom indices
                            self.list_chalcogen.append([self.index_atom_SSE, self.atom_X.index])
                            # append list of distances
                            self.list_chalcogen_distance.append(self.distance)
                            #append list of angle
                            self.list_chalcogen_angle.append([self.angle_theta, self.angle_phi, self.angle_CSX, self.angle_dihedral])
                        
                        # Check hydrogen bond if no hydrogens
                        elif self.MAX_angle_dihedral_chalcogen < self.angle_dihedral and self.angle_CSX <= self.MAX_angle_csx_hbond:
                           
                            # append list of atom indices
                            self.list_hbond.append([self.index_atom_SSE, self.atom_X.index])
                            # append list of distances
                            self.list_hbond_distance.append(self.distance)
                            #append list of angle
                            self.list_hbond_angle.append([self.angle_theta, self.angle_phi, self.angle_CSX, self.angle_dihedral])



            
    #===== Functions =====
    #-----
    def convert_angle_180(self, angle):
        """
        Take angle in degree, in the range 0-360, and convert it to an equavalent range 0-180.
        """
        if angle <= 180:
            return angle
        else:
            return 180 - (angle - 180)

    
    
    #===== Return properties =====    
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if the chalcogen or H-bond, mediated by S/Se, exist between the two residues.
        
        RETURN
            True, 'chalcogen and h-bond'    interaction exist and their is both chalcogen and H-bond
                                            (happend only if the two residues contain S/Se).
            True, 'chalcogen'               interaction exist and is chalcogen.
            True, 'h-bond'                  interaction exists and is H-bond.
            False, False                    no interaction.
        """
        if 0 < len(self.list_chalcogen) and 0 < len(self.list_hbond):
            return True, 'chalcogen and h-bond'
        
        elif 0 < len(self.list_chalcogen):
            return True, 'chalcogen'
        
        elif 0 < len(self.list_hbond):
            return True, 'h-bond'
        
        else:
            return False, False
        
    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return list of atoms index involved in chalogen bond
                       and another list of atoms index involved in H-bond
        RETURN         2 List of list: [[index_S/Se,index_N/O],...] , [[index_S/Se,index_N/O],...]
        """
        return self.list_chalcogen, self.list_hbond
    
    
    @property
    def get_angle(self):
        """
        DESCRIPTION
            Return a list of angles for chalcogen bonds and another for H-bonds
            
        RETURN
            2 List: [[angle_theta, angle_phi, angle_CSX, angle_dihedral],...], [[angle_theta, angle_phi, angle_CSX, angle_dihedral],...]
            
        UNIT
            degree
        """
        return self.list_chalcogen_angle, self.list_hbond_angle
        
    
    @property    
    def get_distance(self):
        """
        DESCRIPTION    Return a list of distance  for chalcogen bonds and another for H-bonds
        RETURN         2 List of list: [distance,...] , [distance,...]
        UNIT           Angstrom
        """
        return self.list_chalcogen_distance, self.list_hbond_distance
        



#=====================================================
#===== Class for S/Se -- Pi interaction
#=====================================================
class sse_aromatic:
    def __init__(self, trajectory, res_index_A, res_index_B, frame=0, MAX_distance=5.5, MIN_pi_angle=60.0, MAX_quadrupole_angle=35.0):
        """
        INTERACTION TYPE    S/Se - aromatic
        SUBTYPE(S)          no

        DESCRIPTION
            Interaction an S, or Se, atom and an aromatic rings.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS        
            MAX_distance          Maximum distance between COM of aromatic ring and the S/Se
                                  Default value: 5.5

            frame                 Frame ID on which to perform the analysis.
                                  Default value: 0: 0
                                  
            MIN_pi_angle            Minimum angle defining the Pi area (the maximum is 90˚).
                                    Default value: 60.0˚
                                    
            MAX_quadrupole_angle    Maximum angle defining the quadrupole area (the maximum is 0˚).
                                    Default value: 35.0˚
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.MAX_distance = MAX_distance
        self.res_index_A = res_index_A
        self.res_index_B = res_index_B
        self.MIN_pi_angle = MIN_pi_angle
        self.MAX_quadrupole_angle = MAX_quadrupole_angle

        
        
        #===== Create dictionaries to  =====
        self.dict_aromatic_ring = {"TYR": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "TRP": 'CD2 CE2 CE3 CZ2 CZ3 CH2',
                                   "PHE": 'CG CD1 CD2 CE1 CE2 CZ',
                                   "HIS": 'CG ND1 CD2 CE1 NE2',
                                   "HID": 'CG ND1 CD2 CE1 NE2',
                                   "HIE": 'CG ND1 CD2 CE1 NE2',
                                   "HSD": 'CG ND1 CD2 CE1 NE2',
                                   "HSE": 'CG ND1 CD2 CE1 NE2',
                                  }
        
        
        self.dict_plane = {"TYR": ['CG',  'CE1',  'CE2'],
                           "TRP": ['CD2', 'CZ2',  'CZ3'],
                           "PHE": ['CG',  'CE1',  'CE2'],
                           "HIS": ['CG',  'CE1',  'NE2'],
                           "HID": ['CG',  'CE1',  'NE2'],
                           "HIE": ['CG',  'CE1',  'NE2'],
                           "HSE": ['CG',  'CE1',  'NE2'],
                           "HSD": ['CG',  'CE1',  'NE2'],
                          }

        
        
        #===== Identify aromatic residue and get it's name =====        
        # Identify the aromatic and the amino AA
        if self.top.residue(self.res_index_A).name in self.dict_aromatic_ring:
            self.res_aromatic = self.res_index_A
        elif self.top.residue(self.res_index_B).name in self.dict_aromatic_ring:
            self.res_aromatic = self.res_index_B
        else:
            raise ValueError('Residues are not TYR TRP PHE HIS')
            
        # Get residues names
        self.res_aromatic_name = self.top.residue(self.res_aromatic).name

        
        #===== Distance between COM of aromatic and S/Se =====
        # Compute COM of aromatic ring
        self.COM_aromatic = md.compute_center_of_mass(self.traj, select=f"resid {self.res_aromatic} and name {self.dict_aromatic_ring[self.res_aromatic_name]}")[0]
        
        # Identify S/Se atom and it's position 
        self.index_atom_SSE = self.top.select(f"resid {self.res_index_A} {self.res_index_B} and element S Se")[0]
        self.coordinate_atom_SSE = self.traj.xyz[0, self.index_atom_SSE]
        
        # create a vector between aromatic COM and SSE & measure it's norm
        self.vector_COM_SSE = Vector.from_points(self.COM_aromatic, self.coordinate_atom_SSE)
        self.distance = self.vector_COM_SSE.norm() *10 # *10 to convert nm to angstrom
        
        
        #===== Create planes and get it's normal vector =====
        # Get position of all atom in aromatic plane
        self.list_aromatic_atom_position = []
        #
        for self.atom in self.dict_plane[self.res_aromatic_name]:
            self.atom_index = self.top.residue(self.res_aromatic).atom(self.atom).index
            self.atom_position = self.traj.xyz[0][self.atom_index]
            self.list_aromatic_atom_position.append(self.atom_position)
        
        # create aromatic plane
        self.plane_aromatic = Plane.from_points(self.list_aromatic_atom_position[0], self.list_aromatic_atom_position[1], self.list_aromatic_atom_position[2])
        
        # Get normal vector od the aromatic plane
        self.vector_normal_plane_aromatic = self.plane_aromatic.normal
        
         
        # ===== Measure angle between the normal vector of aromatic plane and the vector COM --> S/Se =====
        # calculate the angle and convert it to degree
        self.angle = np.rad2deg( self.vector_normal_plane_aromatic.angle_between(self.vector_COM_SSE) )
        
        # Convert the angle to a range of values from 0 to 90 degrees
        if self.angle <= 90:
            None # the angle is already in the range 0 to 90 degrees
        elif self.angle <= 180:
            self.angle = 180 - self.angle
        elif self.angle <= 270:
            self.angle = self.angle - 180
        else:
            self.angle = 360 - self.angle
        
        # The angle is before is between the normal vector of the plan and the vector COM --> S/Se
        # But the normal vector is perpendicular to the plan and geometric criteria are for the angle between the plan and the vector COM --> S/Se.
        # So angle is converted to correspond to be the one betwen the plan and the vector COM --> S/Se.
        self.angle = 90 - self.angle
    
    
    #===== Return results =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if a S/Se interact with an aromatic ring.
        
        RETURN  
            When True, 'S/Se' is replaced by cation or anion.
            True, 'S/Se'-pi              The interaction exist.
            True, 'S/Se'-intermediate    The interaction exist.
            True, 'S/Se'-quadrupole      The interaction exist.
            
            False, False                   The interaction don't exist.
        """
        # get S/Se atom symbol
        symbol = self.top.atom(self.index_atom_SSE).element.symbol
        
        # Take the absolute value of the angle to ensure negative and positive values are considered as the same (ex: -80 is 80)
        if self.distance <= self.MAX_distance and self.MIN_pi_angle <= self.angle <= 90.0:
            return True, f"{symbol}-pi"
        elif self.distance <= self.MAX_distance and self.MAX_quadrupole_angle < self.angle < self.MIN_pi_angle:
            return True, f"{symbol}-intermediate"
        elif self.distance <= self.MAX_distance and 0 <= self.angle <= self.MAX_quadrupole_angle:
            return True, f"{symbol}-quadrupole"
        else:
            return False, False
        
    @property
    def get_angle(self):
        """
        DESCRIPTION    Angle between the S/Se-COM and the normal of the amino plane
        RETURN         angle_planes
        UNIT           degree
        """
        return abs(self.angle)

    @property
    def get_distance(self):
        """
        DESCRIPTION    Distance between COM of aromatic ring and S/Se of the amino group
        RETURN         distance
        UNIT           Angstrom
        """
        return self.distance

    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         
        """
        
        atoms_aromatic = list(self.top.select(f"resid {self.res_aromatic} and name {self.dict_aromatic_ring[self.res_aromatic_name]}"))
        atoms_aromatic.append(self.index_atom_SSE)
        return [atoms_aromatic]




#=====================================================
#===== Function to identify all NB interaction in pair of residues
#=====================================================
def identify_all_interaction_pair(trajectory, pair, frame=0):
    """
    DESCRIPTION
        A fucntion to identify all interaction between a pair of residue.

        It run all interaction commands (with default parameters) and return a
        Pandas DataFrame with informations of each residue, the number of interaction
        doing by the pair, a detail of the interaction type and subtype.
        If the interaction exist: 1, else 0. If the interaction can't exist: NaN.
        
        Exception for salt bridges, it return 1 if bothe h-bond and distance between
        charges are identified, or return 0.5 if only one of them is identified.

        It also return intra-residue C5-Hbond interaction, for each residue.
    
    ARGUMENTS
        pair          pair of residues. The format is a list of MDTraj indices.
                      example: [0,3]

        trajectory    MDTraj trajectory.

    OPTIONAL ARGUMENTS
        frame         trajectory frame.
                      Default value:0
    """
    #===== Create a list to store =====
    list_interaction_pair = []        
        
        
    #----- check interaction: C-bond -----
    try:
        cbond = C_bond(trajectory, pair[0], pair[1], frame=frame)
        
        if cbond.check_interaction == True:
            is_c_bond = 1
        else:
            is_c_bond = 0
    
    except:
        cbond = None
        is_c_bond = np.nan
    
    # add the interaction result to the list_interaction_pair
    list_interaction_pair.append(is_c_bond)
    
    
    
    #----- check interaction: Amino-Pi -----
    try:
        aminoPi = amino_pi(trajectory, pair[0], pair[1], frame=frame)
        
        if aminoPi.check_interaction == True:
            is_amino_pi = 1
        else:
            is_amino_pi = 0
    
    except:
        aminoPi = None
        is_amino_pi = np.nan
    
    # add the interaction result to the list_interaction_pair
    list_interaction_pair.append(is_amino_pi)
    

    
    #----- check interaction: 'Arg-Aromatic' & 'Arg-Arg' -----
    try:
        # get the interaction
        argInvolved = arg_involved(trajectory, pair[0], pair[1], frame=frame)
        argInvolved_type = argInvolved.check_interaction[1]
        argInvolved_subtype = argInvolved.check_interaction[2]
        
        # initialise all variable
        is_arg_aromatic_parallel      = 0
        is_arg_aromatic_perpendicular = 0
        is_arg_aromatic_intermediate  = 0
        is_arg_arg_parallel           = 0
        is_arg_arg_perpendicular      = 0
        is_arg_arg_intermediate       = 0
        
        # check which type of Arg involved interaction is it
        if argInvolved_type == 'Arg-Aromatic' and argInvolved_subtype == 'parallel':
            is_arg_aromatic_parallel = 1
            
        elif argInvolved_type == 'Arg-Aromatic' and argInvolved_subtype == 'perpendicular':
            is_arg_aromatic_perpendicular = 1
            
        elif argInvolved_type == 'Arg-Aromatic' and argInvolved_subtype == 'intermediate':
            is_arg_aromatic_intermediate = 1 
            
        if argInvolved_type == 'Arg-Arg' and argInvolved_subtype == 'parallel':
            is_arg_arg_parallel = 1
            
        elif argInvolved_type == 'Arg-Arg' and argInvolved_subtype == 'perpendicular':
            is_arg_arg_perpendicular = 1
            
        elif argInvolved_type == 'Arg-Arg' and argInvolved_subtype == 'intermediate':
            is_arg_arg_intermediate = 1 
    
    except:
        argInvolved = None
        is_arg_aromatic_parallel      = np.nan
        is_arg_aromatic_perpendicular = np.nan
        is_arg_aromatic_intermediate  = np.nan
        is_arg_arg_parallel           = np.nan
        is_arg_arg_perpendicular      = np.nan
        is_arg_arg_intermediate       = np.nan
    
    # add the interaction result to the list_interaction_pair
    list_interaction_pair.extend([is_arg_aromatic_parallel, is_arg_aromatic_perpendicular, is_arg_aromatic_intermediate, \
                                 is_arg_arg_parallel, is_arg_arg_perpendicular, is_arg_arg_intermediate])
    
    
    
    #----- check interaction: Aromatic-Aromatic -----
    try:
        # get the interaction
        aromaticAromatic = aromatic_aromatic (trajectory, pair[0], pair[1], frame=frame)
        
        #initialise result variable
        is_aromatic_aromatic_parallel = 0
        is_aromatic_aromatic_offset   = 0
        is_aromatic_aromatic_coplanar = 0
        is_aromatic_aromatic_Yshaped  = 0
        is_aromatic_aromatic_Tshaped  = 0
        is_aromatic_aromatic_intermediate = 0
        
        # check which type of aromatic-aromatic interaction is it
        if aromaticAromatic.check_interaction[1] == 'parallel':
            is_aromatic_aromatic_parallel = 1
            
        elif aromaticAromatic.check_interaction[1] == 'offset':
            is_aromatic_aromatic_offset = 1
            
        elif aromaticAromatic.check_interaction[1] == 'coplanar':
            is_aromatic_aromatic_coplanar = 1
            
        elif aromaticAromatic.check_interaction[1] == 'Y-shaped':
            is_aromatic_aromatic_Yshaped = 1
            
        elif aromaticAromatic.check_interaction[1] == 'T-shaped':
            is_aromatic_aromatic_Tshaped = 1

        elif aromaticAromatic.check_interaction[1] == 'intermediate':
            is_aromatic_aromatic_intermediate = 1
    
    except:
        aromaticAromatic = None
        is_aromatic_aromatic_parallel = np.nan
        is_aromatic_aromatic_offset   = np.nan
        is_aromatic_aromatic_coplanar = np.nan
        is_aromatic_aromatic_Yshaped  = np.nan
        is_aromatic_aromatic_Tshaped  = np.nan
        is_aromatic_aromatic_intermediate = np.nan
    
    # add the interaction result to the list_interaction_pair
    list_interaction_pair.extend([is_aromatic_aromatic_parallel, is_aromatic_aromatic_offset, is_aromatic_aromatic_coplanar, \
                                 is_aromatic_aromatic_Yshaped, is_aromatic_aromatic_Tshaped, is_aromatic_aromatic_intermediate])
    
    
    
    #----- check interaction: charge-aromatic -----
    try:
        # get the interaction
        chargeAromatic = charge_aromatic (trajectory, pair[0], pair[1], frame=frame)

        #initialise result variable
        is_cation_pi = 0
        is_cation_intermediate = 0
        is_cation_quadrupole = 0
        is_anion_pi = 0
        is_anion_intermediate = 0
        is_anion_quadrupole = 0

        # check which type is it
        if chargeAromatic.check_interaction[1] == 'cation-pi':
            is_cation_pi = 1

        elif chargeAromatic.check_interaction[1] == 'cation-intermediate':
            is_cation_intermediate = 1

        elif chargeAromatic.check_interaction[1] == 'cation-quadrupole':
            is_cation_quadrupole = 1

        elif chargeAromatic.check_interaction[1] == 'anion-pi':
            is_anion_pi = 1

        elif chargeAromatic.check_interaction[1] == 'anion-intermediate':
            is_anion_intermediate = 1

        elif chargeAromatic.check_interaction[1] == 'anion-quadrupole':
            is_anion_quadrupole = 1
    
    except:
        chargeAromatic = None
        is_cation_pi = np.nan
        is_cation_intermediate = np.nan
        is_cation_quadrupole = np.nan
        is_anion_pi = np.nan
        is_anion_intermediate = np.nan
        is_anion_quadrupole = np.nan
    
    # add the interaction result to the list_interaction_pair
    list_interaction_pair.extend([is_cation_pi, is_cation_intermediate, is_cation_quadrupole, \
                                 is_anion_pi, is_anion_intermediate, is_anion_quadrupole])
    
    
    
    #----- check interaction: Clash or repulsion between charges -----
    try:
        # get the interaction
        chargeClashRepulsion = charge_clash_repulsion (trajectory, pair[0], pair[1], frame=frame)
        
        # initalise varaibale
        is_charge_repulsion = 0
        is_charge_clash = 0
        
        # check which type is it
        if chargeClashRepulsion.check_interaction[2] == 'repulsion':
            is_charge_repulsion = 1

        elif chargeClashRepulsion.check_interaction[2] == 'clash':
            is_charge_clash = 1
    
    except:
        chargeClashRepulsion = None
        is_charge_repulsion = np.nan
        is_charge_clash = np.nan
    
    # add the interaction result to the list_interaction_pair
    list_interaction_pair.extend([is_charge_repulsion, is_charge_clash])
    
    
    
    #----- check interaction: hydrogen bonds -----
    try:
        hbond = hydrogen_bond(trajectory, pair[0], pair[1], frame=frame)
        
        if hbond.check_interaction == True:
            is_hydrogen_bond = 1
        else:
            is_hydrogen_bond = 0
    
    except:
        hbond = None
        is_hydrogen_bond = np.nan
    
    # add the interaction result to the list_interaction_pair
    list_interaction_pair.append(is_hydrogen_bond)
    
    
    
    #----- check interaction: hydrophobic interaction & hydrophobe/hydrophile clash -----
    try:
        # get the interaction
        hydrophobic_ClashRepulsion = hydrophobic(trajectory, pair[0], pair[1], frame=frame)
        
        # initialise variable
        is_hydrophobic_interaction = 0
        is_hydrophobe_hydrophile_clash = 0
        is_hydrophobe_hydrophile_repulsion = 0
        
        # check which type is it
        if hydrophobic_ClashRepulsion.check_interaction[2] == 'hydrophobic':
            is_hydrophobic_interaction = 1
        
        elif hydrophobic_ClashRepulsion.check_interaction[2] == 'clash':
            is_hydrophobe_hydrophile_clash = 1

        elif hydrophobic_ClashRepulsion.check_interaction[2] == 'repulsion':
            is_hydrophobe_hydrophile_repulsion = 1
    
    except:
        hydrophobic_ClashRepulsion = None
        is_hydrophobic_interaction = np.nan
        is_hydrophobe_hydrophile_clash = np.nan
        is_hydrophobe_hydrophile_repulsion = np.nan


    # add the interaction result to the list_interaction_pair
    list_interaction_pair.extend([is_hydrophobic_interaction, is_hydrophobe_hydrophile_clash, is_hydrophobe_hydrophile_repulsion])
    
    
    
    #----- check interaction: n-->pi* -----
    try:
        # get the interaction
        npi = n_pi(trajectory, pair[0], pair[1], frame=frame)
        
        # initialise result variable
        is_n_pi_regular = 0
        is_n_pi_reciprocal = 0
        
        # search which type is it
        if npi.check_interaction[1] == 'regular':
            is_n_pi_regular = 1
        
        elif npi.check_interaction[1] == 'reciprocal':
            is_n_pi_reciprocal = 1  
    
    except:
        npi = None
        is_n_pi_regular = np.nan
        is_n_pi_reciprocal = np.nan
    
    # add the interaction result to the list_interaction_pair
    list_interaction_pair.extend([is_n_pi_regular, is_n_pi_reciprocal])
    
    
    
    #----- check interaction: pi H-bond -----
    try:
        piHbond = pi_hbond(trajectory, pair[0], pair[1], frame=frame)
        
        if piHbond.check_interaction == True:
            is_pi_hbond = 1
        else:
            is_pi_hbond = 0

    except:
        piHbond = None
        is_pi_hbond = np.nan

    
    # add the interaction result to the list_interaction_pair
    list_interaction_pair.append(is_pi_hbond)  
    
    
    
    #----- check interaction: Salt bridge -----
    try:
        # get the interaction
        saltbridge = salt_bridge(trajectory, pair[0], pair[1], frame=frame)
        
        # Salt bridge is 1 if both charge distance & H-bond
        if saltbridge.check_interaction[0] == True and saltbridge.check_interaction[1] == True:
            is_salt_bridge = 1
            
        # Salt bridge is 0.5 if only charge distance or only H-bond  
        elif saltbridge.check_interaction[0] == True or saltbridge.check_interaction[1] == True:
            is_salt_bridge = 0.5
            
        else:
            is_salt_bridge = 0

    except:
        saltbridge = None
        is_salt_bridge = np.nan

    # add the interaction result to the list_interaction_pair
    list_interaction_pair.append(is_salt_bridge)
    
    
    
    #----- check interaction: S/Se-Aromatic -----
    try:
        # get the interaction
        SSeAromatic = sse_aromatic(trajectory, pair[0], pair[1], frame=frame)

        # initalise variable
        is_S_pi = 0
        is_S_intermediate = 0
        is_S_quadrupole = 0
        is_Se_pi = 0
        is_Se_intermediate = 0
        is_Se_quadrupole = 0
        
        # search which type is it
        if SSeAromatic.check_interaction[1] == 'S-pi':
            is_S_pi = 1
            
        elif SSeAromatic.check_interaction[1] == 'S-intermediate':
            is_S_intermediate = 1
            
        elif SSeAromatic.check_interaction[1] == 'S-quadrupole':
            is_S_quadrupole = 1
            
        elif SSeAromatic.check_interaction[1] == 'Se-pi':
            is_Se_pi = 1
            
        elif SSeAromatic.check_interaction[1] == 'Se-intermediate':
            is_Se_intermediate = 1
            
        elif SSeAromatic.check_interaction[1] == 'Se-quadrupole':
            is_Se_quadrupole = 1            

    except:
        SSeAromatic = None
        is_S_pi = np.nan
        is_S_intermediate = np.nan
        is_S_quadrupole = np.nan
        is_Se_pi = np.nan
        is_Se_intermediate = np.nan
        is_Se_quadrupole = np.nan

    # add the interaction result to the list_interaction_pair
    list_interaction_pair.extend([is_S_pi, is_S_intermediate, is_S_quadrupole, \
                                 is_Se_pi, is_Se_intermediate, is_Se_quadrupole])
        
    
    
    #----- check interaction: S/Se mediated H-bon & chalcogen -----
    try:
        # get the interaction
        SSE_hbond_chalcogen = sse_hydrogen_chalcogen_bond(trajectory, pair[0], pair[1], frame=frame)

        # inialize variable
        is_sse_hbond = 0
        is_sse_chalcogen = 0

        if SSE_hbond_chalcogen.check_interaction[1] == 'h-bond':
            is_sse_hbond = 1
        elif SSE_hbond_chalcogen.check_interaction[1] == 'chalcogen':
            is_sse_chalcogen = 1

    except:
        SSE_hbond_chalcogen = None
        is_sse_hbond = np.nan
        is_sse_chalcogen = np.nan

    # add the interaction result to the list_interaction_pair
    list_interaction_pair.extend([is_sse_hbond, is_sse_chalcogen])
    
    
    
    #----- check interaction: van der Waals -----
    try:
        vdw = van_der_waals(trajectory, pair[0], pair[1], frame=frame)
        # get number of vdw contacts
        van_der_waals_nb_contacts = vdw.get_number_contacts

    except:
        vdw = None
        van_der_waals_nb_contacts = np.nan
        is_van_der_waals = np.nan

    if 0 < van_der_waals_nb_contacts:
        is_van_der_waals = 1
    else:
        is_van_der_waals = 0

    # add the interaction result to the list_interaction_pair
    list_interaction_pair.append(is_van_der_waals)


    #===== Perform intra-residue C5 Hbonds =====
    #----- check if residues 1 and 2 have C5-Hbond (intra-residue H-bond) -----
        #..... Residue 1
    try:
        # search C5-Hbond
        residue_1_C5Hbond = C5_hydrogen_bond(trajectory, pair[0], frame=frame)
        
        # check if exist for residue 1 or not
        if residue_1_C5Hbond.check_interaction == True:
            is_residue_1_C5Hbond = 1
        else:
            is_residue_1_C5Hbond = 0
    
    except:
        residue_1_C5Hbond = None
        is_residue_1_C5Hbond = np.nan
    
    #..... Residue 2
    try:
        # search C5-Hbond
        residue_2_C5Hbond = C5_hydrogen_bond(trajectory, pair[1], frame=frame)
        
        # check if exist for residue 2 or not
        if residue_2_C5Hbond.check_interaction == True:
            is_residue_2_C5Hbond = 1
        else:
            is_residue_2_C5Hbond = 0
    
    except:
        residue_2_C5Hbond = None
        is_residue_2_C5Hbond = np.nan 

    
    #===== return result =====
    #----- Get information in the topology, concerning the two residues -----
    # for residue 1
    residue_1 = trajectory.topology.residue(pair[0])
    residue_1_resseq = residue_1.resSeq
    residue_1_name   = residue_1.name
    residue_1_chain  = residue_1.chain.index
    
    # for residue 2
    residue_2 = trajectory.topology.residue(pair[1])
    residue_2_resseq = residue_2.resSeq
    residue_2_name   = residue_2.name
    residue_2_chain  = residue_2.chain.index
    
    # check if the two residues are consecutive or not: residue i and residue i+1
    residues_index_difference = abs(pair[0] - pair[1])
    
    if residue_1_chain == residue_2_chain and residues_index_difference == 1:
        consecutive = 1
    else:
        consecutive = 0
    
    
    #----- calculate the number of interaction in the pair -----
    #Check if the list_interaction_pair contain at least one 1 or one 0.5
    array_interaction_pair = np.array(list_interaction_pair)

    # np.ceil is used to up 0.5 to 1 for the salt bridge
    # np.nansum is used to calculate the sum without taking in account the np.nan
    number_interactions = np.nansum( np.ceil( array_interaction_pair ) )
    
    #----- create a pandas dataframe of interactions -----
    df = pd.DataFrame([{"residue_1_chain" : residue_1_chain,
                        "residue_1_index" : pair[0],
                        "residue_1_resSeq": residue_1_resseq,
                        "residue_1_name"  : residue_1_name,
                        #
                        "residue_2_chain" : residue_2_chain,
                        "residue_2_index" : pair[1],
                        "residue_2_resSeq": residue_2_resseq,
                        "residue_2_name"  : residue_2_name,
                        #
                        "consecutive_residues": consecutive,
                        "number_interactions" : number_interactions,
                        #
                        "residue_1_C5Hbond": is_residue_1_C5Hbond,
                        "residue_2_C5Hbond": is_residue_2_C5Hbond,
                        #
                        "c_bond" : is_c_bond,
                        #
                        "amino_pi" : is_amino_pi,
                        #
                        "arg_aromatic_parallel"      : is_arg_aromatic_parallel,    
                        "arg_aromatic_perpendicular" : is_arg_aromatic_perpendicular,
                        "arg_aromatic_intermediate " : is_arg_aromatic_intermediate,
                        "arg_arg_parallel"      : is_arg_arg_parallel,     
                        "arg_arg_perpendicular" : is_arg_arg_perpendicular,     
                        "arg_arg_intermediate"  : is_arg_arg_intermediate, 
                        #
                        "aromatic_aromatic_parallel" : is_aromatic_aromatic_parallel,
                        "aromatic_aromatic_offset"   : is_aromatic_aromatic_offset,  
                        "aromatic_aromatic_coplanar" : is_aromatic_aromatic_coplanar,
                        "aromatic_aromatic_Yshaped"  : is_aromatic_aromatic_Yshaped, 
                        "aromatic_aromatic_Tshaped"  : is_aromatic_aromatic_Tshaped,
                        "aromatic_aromatic_intermediate" : is_aromatic_aromatic_intermediate, 
                        #
                        "cation_pi"           : is_cation_pi,
                        "cation_intermediate" : is_cation_intermediate,
                        "cation_quadrupole"   : is_cation_quadrupole,
                        "anion_pi"            : is_anion_pi,
                        "anion_intermediate"  : is_anion_intermediate,
                        "anion_quadrupole"    : is_anion_quadrupole,
                        #
                        "charge_repulsion" : is_charge_repulsion,
                        "charge_clash"     : is_charge_clash,
                        #
                        "hydrogen_bond" : is_hydrogen_bond,
                        #
                        "hydrophobic_interaction"     : is_hydrophobic_interaction,
                        "hydrophobe_hydrophile_clash" : is_hydrophobe_hydrophile_clash,
                        "hydrophobe_hydrophile_repulsion" : is_hydrophobe_hydrophile_repulsion,
                        #
                        "n_pi_regular"    : is_n_pi_regular,
                        "n_pi_reciprocal" : is_n_pi_reciprocal,
                        #
                        "pi_hbond" : is_pi_hbond,
                        #
                        "salt_bridge" : is_salt_bridge,
                        #
                        "S_pi"            : is_S_pi,
                        "S_intermediate"  : is_S_intermediate,
                        "S_quadrupole"    : is_S_quadrupole,
                        "Se_pi"           : is_Se_pi,
                        "Se_intermediate" : is_Se_intermediate,
                        "Se_quadrupole"   : is_Se_quadrupole,
                        #
                        "sse_hbond"     : is_sse_hbond,
                        "sse_chalcogen" : is_sse_chalcogen,
                        #
                        "van_der_waals_nb_contacts" : van_der_waals_nb_contacts,
    }])
    
    #----- create a pandas dataframe of resulting classes (python object) -----
    df_class = pd.DataFrame({ 
        "residue_1_index" : pair[0],
        "residue_2_index" : pair[1],
        "residue_1_C5Hbond" : residue_1_C5Hbond,
        "residue_2_C5Hbond" : residue_2_C5Hbond,
        "c_bond"   : cbond,
        "amino_pi" : aminoPi,
        "arg_involved" : argInvolved,
        "aromatic_aromatic" : aromaticAromatic,
        "charge_aromatic"   : chargeAromatic,
        "charge_clash_repulsion" : chargeClashRepulsion,
        "hydrogen_bond" : hbond,
        "hydrophobic"   : hydrophobic_ClashRepulsion,
        "n_pi" : npi,
        "pi_hydrogen_bond" : piHbond,
        "salt_bridge"  : saltbridge,
        "sse_aromatic" : SSeAromatic,
        "sse_hydrogen_chalcogen_bond" : SSE_hbond_chalcogen,
        "van_der_waals" : vdw,
    }, index=[0])

    return df, df_class





#=====================================================
#===== Function to identify all NB interaction in a complete structure
#=====================================================

#----- Function used by multiprocessing -----
def task(trajectory, pair, frame, path_class_outfile, write_class_outfile):
    """
    Is function is defined to be used by 'interaction_table_whole_system()' with multiprocessing.
    Not designed to be used by user !
    """
    #----- Check which interaction types are performed between the two residues -----
    interactions = identify_all_interaction_pair(trajectory, pair, frame=frame)
    df_interactions = interactions[0]
    df_class_interactions = interactions[1]

    # if their is not interaction, pass an go to next item in 'iterable'
    if not isinstance(df_interactions, pd.DataFrame) or not isinstance(df_class_interactions, pd.DataFrame):
        return None, None

    else:
        #----- Calculate Coulomb and LJ energies using AMBER and CHARMM -----
        try:
            energy = omm_coulomb_lj(trajectory, pair[0], pair[1], frame=frame)
            #
            energy_total_amber   = energy.get_energy[0]
            energy_coulomb_amber = energy.get_energy_coulomb[0]
            energy_lj_amber      = energy.get_energy_LJ[0]
            #
            energy_total_charmm   = energy.get_energy[1]
            energy_coulomb_charmm = energy.get_energy_coulomb[1]
            energy_lj_charmm      = energy.get_energy_LJ[1]

        except:
            energy = None
            #
            energy_total_amber   = np.nan
            energy_coulomb_amber = np.nan
            energy_lj_amber      = np.nan
            #
            energy_total_charmm   = np.nan
            energy_coulomb_charmm = np.nan
            energy_lj_charmm      = np.nan   

            
        #----- Add energies and C5-Hbond to the 'df_interactions' row -----
        df_interactions["energy_total_charmm"] = energy_total_charmm
        df_interactions["energy_coulomb_charmm"] = energy_coulomb_charmm
        df_interactions["energy_lj_charmm"] = energy_lj_charmm
        df_interactions["energy_total_amber"] = energy_total_amber
        df_interactions["energy_coulomb_amber"] = energy_coulomb_amber
        df_interactions["energy_lj_amber"] = energy_lj_amber

        if write_class_outfile == True:
            class_outfile = f"{path_class_outfile}_{pair[0]}_{pair[1]}.pkl.gz"
            df_class_interactions.to_pickle(class_outfile)

        return df_interactions
        


#----- Function to analyse the whole system -----
def interaction_table_whole_system(trajectory, list_pairs="all", frame=0, MAX_CA_distance=14.0, use_tqdm=False, \
                                   write_outfile=True, path_table_outfile="interaction_table_whole_system.csv", \
                                   write_class_outfile=False, path_class_outfile="class_"):
    """
    DESCRIPTION
        In a system (protein/complex) analyse all possible pairs (or given pairs) and return a complete table of their interaction types,
        with: 
            - Coulomb and Lennard-Jones energyes as calculated by openMM using both AMBER and CHARMM,

    ARGUMENTS
        trajectory    MDTraj trajectory.

    OPTIONAL ARGUMENTS
        frame               Frame ID on which to perform the analysis.
                            Default value: 0
        
        list_pairs          List of pair of residue. By default take all possible pair in the system.
                            Else can be set by user, using this format: [[id1,id2], [id3,id4], ...]
        
        MAX_CA_distance     Maximum distance between CA atom of the two residue.
                            Pairs with a greater distance will be ignored.
                            Default value: 14.0 angstrom

        use_tqdm            Display tqdm proggress bar (True/False).
                            Default value: False

        write_outfile       Write the complete interaction table in CSV format.
                            (True/False)
                            Default value: True

        path_table_outfile   Path, containing the name, of the exported CSV table.
                             By default the file is exported in the curent location
                             and the name: "interaction_table_whole_system.csv"
                             Custom: "my/path/interaction_table_whole_system.csv"

        write_class_outfile  Write all computed interaction classes in pickle format,
                             compressed in gz.
                             (True/False)
                             Default value: False

        path_class_outfile   Path, without the name, of the exported computed classes.
                             By default the files are exported in the curent location
                             not the name, nor the file format.
                             Custom: "my/path/system_class"
    """
    #===== Initialize trajectory variable =====
    trajectory = trajectory[frame]
    
    
    #===== Get the list of ressidue pair =====
    if list_pairs == "all":
        list_resi = [i.index for i in trajectory.topology.residues]
        list_residue_pairs = [tuple(sorted(i)) for i in list(combinations(list_resi , 2))]

    else:
        list_residue_pairs = list_pairs

        
    #===== Select only pairs with CA-CA distance <= MAX_CA_distance =====
    # use MDTraj to compute distances between all pairs
    all_distances_indices = md.compute_contacts(trajectory, contacts=list_residue_pairs, scheme='ca')
    
    # extract distance and pairs ID
    distances = all_distances_indices[0][0] *10   #*10 to convert nm to angstrom
    pairs     = all_distances_indices[1]

    # Create a boolean mask for distances <= MAX_CA_distance
    mask = distances <= MAX_CA_distance

    # Use the mask to filter pairs arrays
    filtered_pairs = pairs[mask]


    #===== create the process pool =====    
    with Pool() as pool:
        
        # prepare arguments
        items = [(trajectory, i, frame, path_class_outfile, write_class_outfile) for i in filtered_pairs]
        
        # loop over all residue pair to analyse their interactions
        if use_tqdm == True:
            results = pool.starmap(task, tqdm(items, total=len(items)))
        else:
            results = pool.starmap(task, items)


    
    #===== Return results =====
    # Initialize final tables
    df = pd.DataFrame() # table containing the analysis result summary
    
    # Append the final tables with the row 
    df = pd.concat([i[0] for i in results if isinstance(i[0], pd.DataFrame)], ignore_index=True)
    
    if write_outfile == True:
        df.to_csv(path_table_outfile, index=False)
    
    # return the complete dataframe
    return df
    




#=====================================================
#===== Function to identify interaction area
#=====================================================
class locate:
    def __init__ (self, result, directory='.', file_name_structure=None):
        """
        DESCRIPTION
            Take an interaction type class object and locate the interaction area where the interaction take place:
            backbone-backbone (BB-BB), backbone-sidechain (name-BB), sidechain-sidechain (name1-name2).

            If a cleaned structure file is provided, the class can return 'code_complete'
            and 'code_name_secondary_structure' for backbone and side chaine.

            Warning: It always return result, also if the interaction is "False", because atom's
                     location/area (sidechain or bakbone) is independant on the interaction existance.
            Warning: Error "AttributeError: 'NoneType' object has no attribute 'get_atoms'" append whent the interaction type
                     class is None, i.e. don't exist.

            
        ARGUMENTS
            result    interaction type class object  
            
        OPTIONAL ARGUMENTS
            directory   Directry whre the CSV file of cleaned structure information is located.
                        Default vale: '.'
                        
            file_name_structure    Full name of motif in the name of the cleaned structure file.
                                   Default vale: None            
        """
        self.result = result
        self.directory = directory
        self.file_name_structure = file_name_structure

        
        #===== Try to find secondary structure and protein region information to make dictionnaries =====
        try:
            # Set file name
            if self.file_name_structure == None:
                self.file_name_structure = 'clean_structure'
            
            # Read as dataframe
            self.file_structure = glob.glob(os.path.join(self.directory, f'*{self.file_name_structure}*.csv'))[0]
            self.df_structure   = pd.read_csv(self.file_structure)
            
            # Create dictionnaries
            self.dict_residue_index2codeComplete = self.df_structure.set_index('index')['code_complete'].to_dict()
            self.dict_residue_index2codeSimple   = self.df_structure.set_index('index')['code_name_secondary_structure'].to_dict()
            self.dict_residue_index2ss           = self.df_structure.set_index('index')['secondary_structure'].to_dict()
            self.dict_residue_index2pr           = self.df_structure.set_index('index')['protein_region'].to_dict()
                
        except:
            self.dict_residue_index2codeComplete = None
            self.dict_residue_index2codeSimple   = None
            self.dict_residue_index2ss           = None
            self.dict_residue_index2pr           = None
        
        
        #===== Get indices of residues =====
        # depending on the class the residue indices definid in the self can be self.res_A/B or self.res_index_A/B
        try:
            self.residue_index_A = self.result.res_A
            self.residue_index_B = self.result.res_B
        except:
            try:
                self.residue_index_A = self.result.res_index_A
                self.residue_index_B = self.result.res_index_B
            except:
                self.residue_index_A = self.result.res_R_index
                self.residue_index_B = self.result.res_aromatic_index

        self.pair_index = '_'.join(map(str, list(sorted([self.residue_index_A, self.residue_index_B]))))
    
        
        #===== Get residues information =====
        # residue names
        self.residue_name_A = self.result.top.residue(self.residue_index_A).name
        self.residue_name_B = self.result.top.residue(self.residue_index_B).name
    
        # code complete
        if self.dict_residue_index2codeComplete != None:
            self.residue_codeComplete_A = self.dict_residue_index2codeComplete[self.residue_index_A]
            self.residue_codeComplete_B = self.dict_residue_index2codeComplete[self.residue_index_B]
        else:
            self.residue_codeComplete_A = None
            self.residue_codeComplete_B = None
    
        # code simple
        if self.dict_residue_index2codeSimple != None:
            self.residue_codeSimple_A = self.dict_residue_index2codeSimple[self.residue_index_A]
            self.residue_codeSimple_B = self.dict_residue_index2codeSimple[self.residue_index_B]
        else:
            self.residue_codeSimple_A = None
            self.residue_codeSimple_B = None
    
        # code backbone
        if self.dict_residue_index2ss != None:
            self.backbone_codeSimple_A   = f"BB_{self.dict_residue_index2ss[self.residue_index_A]}"
            self.backbone_codeSimple_B   = f"BB_{self.dict_residue_index2ss[self.residue_index_B]}"
    
            if self.dict_residue_index2pr != None:
                self.backbone_codeComplete_A = f"BB_{self.dict_residue_index2ss[self.residue_index_A]}_{self.dict_residue_index2pr[self.residue_index_A]}"
                self.backbone_codeComplete_B = f"BB_{self.dict_residue_index2ss[self.residue_index_B]}_{self.dict_residue_index2pr[self.residue_index_B]}"
            else:
                self.backbone_codeComplete_A = None
                self.backbone_codeComplete_B = None
        
        else:
            self.backbone_codeSimple_A   = None
            self.backbone_codeSimple_B   = None
            self.backbone_codeComplete_A = None
            self.backbone_codeComplete_B = None
        
        
        #===== Get set of atoms in the sidechain and backbone of residues =====
        # for residue A
        self.set_sidechain_residue_A = set(self.result.top.select(f"sidechain and resid {self.residue_index_A}"))
        self.set_backone_residue_A   = set(self.result.top.select(f"backbone and resid {self.residue_index_A}"))
        self.set_atoms_residue_A     = set(self.result.top.select(f"resid {self.residue_index_A}"))
        
        # for residue B
        self.set_sidechain_residue_B = set(self.result.top.select(f"sidechain and resid {self.residue_index_B}"))
        self.set_backone_residue_B   = set(self.result.top.select(f"backbone and resid {self.residue_index_B}"))
        self.set_atoms_residue_B     = set(self.result.top.select(f"resid {self.residue_index_B}"))
        
        # Define set for all sidechains and backbones 
        self.set_atoms_sidechain = self.set_sidechain_residue_A.union(self.set_sidechain_residue_B) 
        self.set_atoms_backbone  = self.set_backone_residue_A.union(self.set_backone_residue_B)
        
        
        #===== Get list of atoms performing the interaction =====
        #----- set list to store results -----
        self.list_interaction_area          = []
        self.list_interaction_code_name     = []
        self.list_interaction_code_complete = []
        self.list_interaction_code_simple   = []
        
        #----- Read the .get_atoms output -----
        # if the output is a tupple of list (occure when: return list1, list2) merge the lists
        if isinstance(self.result.get_atoms, tuple):
            self.atoms_list = []
            for x in self.result.get_atoms:
                self.atoms_list.extend(x)

        # if the output is a list only take it as atom list
        else:
            self.atoms_list = self.result.get_atoms



        # correction for pi_hbond interaction
        if self.result.__class__.__name__.lower() == 'pi_hbond':
            newlist = []
            
            for H, X in self.atoms_list:
                if H in self.set_atoms_residue_A:
                    newlist.append(self.set_sidechain_residue_B.union([H]))
                
                elif H in self.set_atoms_residue_B:
                    newlist.append(self.set_sidechain_residue_A.union([H]))
            
            self.atoms_list = newlist



        for self.atom_indices in self.atoms_list:
            self.set_atom_interaction = set(self.atom_indices)
            
            #----- Search for protein area -----
            if any(self.set_atom_interaction & self.set_atoms_sidechain) and not any(self.set_atom_interaction & self.set_atoms_backbone):
                # Add area
                self.list_interaction_area.append("sidechain-sidechain")
    
                # Add code name
                self.code_name = '-'.join(sorted([self.residue_name_A, self.residue_name_B]))
                self.list_interaction_code_name.append(self.code_name)
    
                # Add code complete
                if self.residue_codeComplete_A and self.residue_codeComplete_B:
                    self.code_complete = '-'.join(sorted([self.residue_codeComplete_A, self.residue_codeComplete_B]))
                    self.list_interaction_code_complete.append(self.code_complete)
                else:
                    self.list_interaction_code_complete.append("NaN-NaN")
                
                # Add code simple
                if self.residue_codeSimple_A and self.residue_codeSimple_B:    
                    self.code_simple = '-'.join(sorted([self.residue_codeSimple_A, self.residue_codeSimple_B]))
                    self.list_interaction_code_simple.append(self.code_simple)
                else:
                    self.list_interaction_code_simple.append("NaN-NaN")
    
            
            #..... check for backbone-backbone .....
            elif not any(self.set_atom_interaction & self.set_atoms_sidechain) and any(self.set_atom_interaction & self.set_atoms_backbone):
                # Add area
                self.list_interaction_area.append("backbone-backbone")
    
                # Add code name
                self.list_interaction_code_name.append("BB-BB")
    
                # Add code complete
                if self.backbone_codeComplete_A and self.backbone_codeComplete_B:
                    self.code_complete = '-'.join(sorted([self.backbone_codeComplete_A, self.backbone_codeComplete_B]))
                    self.list_interaction_code_complete.append(self.code_complete)
                else:
                    self.list_interaction_code_complete.append("NaN-NaN")
                
                # Add code simple
                if self.backbone_codeSimple_A and self.backbone_codeSimple_B:    
                    self.code_simple = '-'.join(sorted([self.backbone_codeSimple_A, self.backbone_codeSimple_B]))
                    self.list_interaction_code_simple.append(self.code_simple)
                else:
                    self.list_interaction_code_simple.append("NaN-NaN")
    
            
            #..... backbone-sidechain .....
            elif any(self.set_atom_interaction & self.set_atoms_sidechain) and any(self.set_atom_interaction & self.set_atoms_backbone):
    
                # Idendify which sidechain interaction with backbone
                if any(self.set_atom_interaction & self.set_sidechain_residue_A) and any(self.set_atom_interaction & self.set_backone_residue_B):
                    # residue A is side chain and ressidue B is backbone
                    self.sidechain_name         = self.residue_name_A
                    self.sidechain_codeComplete = self.residue_codeComplete_A
                    self.sidechain_codeSimple   = self.residue_codeSimple_A
                    self.backbone_codeSimple    = self.backbone_codeSimple_B
                    self.backbone_codeComplete  = self.backbone_codeComplete_B
    
                elif any(self.set_atom_interaction & self.set_sidechain_residue_B) and any(self.set_atom_interaction & self.set_backone_residue_A):
                    # residue B is side chain and ressidue A is backbone
                    self.sidechain_name         = self.residue_name_B
                    self.sidechain_codeComplete = self.residue_codeComplete_B
                    self.sidechain_codeSimple   = self.residue_codeSimple_B
                    self.backbone_codeSimple    = self.backbone_codeSimple_A
                    self.backbone_codeComplete  = self.backbone_codeComplete_A
                
                # Add area
                self.list_interaction_area.append("backbone-sidechain")
    
                # Add code name
                self.list_interaction_code_name.append(f"{self.sidechain_name}-BB")
    
                # Add code complete
                if self.sidechain_codeComplete and self.backbone_codeComplete:
                    self.list_interaction_code_complete.append(f"{self.sidechain_codeComplete}-{self.backbone_codeComplete}")
                else:
                    self.list_interaction_code_complete.append("NaN-NaN")
                
                # Add code simple
                if self.sidechain_codeSimple and self.backbone_codeSimple:
                    self.list_interaction_code_simple.append(f"{self.sidechain_codeSimple}-{self.backbone_codeSimple}")
                else:
                    self.list_interaction_code_simple.append("NaN-NaN")
    
    
    #===== Return informations =====
    @property
    def get_interaction_type(self):
        """
        DESCRIPTION    Return interaction type name. (It's the class name)
        RETURN         result.__class__.__name__.lower()
        """
        return self.result.__class__.__name__.lower()
        
    @property
    def get_atoms(self):
        """
        DESCRIPTION    Return atoms invoved in the interaction.
        RETURN         result.atoms_list
        """
        return self.atoms_list
        
    @property
    def get_area(self):
        """
        DESCRIPTION    Return area list corresponding to atoms item in '.get_atoms' list.
        RETURN         list_interaction_area
        """
        return self.list_interaction_area
    
    @property
    def get_code(self):
        """
        DESCRIPTION    Return interaction code name, code name secondary structure (= code simple), and
                       code complete corresponding to atoms item in '.get_atoms' list.
        RETURN         list_interaction_code_name, list_interaction_code_simple, list_interaction_code_complete
        """
        return self.list_interaction_code_name, self.list_interaction_code_simple, self.list_interaction_code_complete


    @property
    def get_table_interaction(self):
        """
        DESCRIPTION    Return the count for interaction type areas.
                       'pair_index' is set as index of the table to easly merge with other interaction type row (or table).
        RETURN         pandas dataframe
        """
        # Sample list of elements
        elements = [self.result.__class__.__name__.lower() + "_" + x for x in self.list_interaction_area]
        
        # Count the occurrences of each element
        element_counts = Counter(elements)
        
        # Convert the Counter object to a Pandas DataFrame
        df = pd.DataFrame(list(element_counts.items()), columns=['pair_index', self.pair_index])
        df.set_index('pair_index', inplace=True)
        
        return df.T


    @property
    def get_table_code_name(self):
        """
        DESCRIPTION    Return the count of code_name pair for interaction type.
                       'code_name' is set as index of the table to easly merge with other interaction type row (or table).
        RETURN         pandas dataframe
        """        
        # Count the occurrences of each element
        element_counts = Counter(self.list_interaction_code_name)
        
        # Convert the Counter object to a Pandas DataFrame
        df = pd.DataFrame(list(element_counts.items()), columns=['pair_code_name', self.result.__class__.__name__.lower()])
        df.set_index('pair_code_name', inplace=True)
        
        return df

    
    
    @property
    def get_table_code_simple(self):
        """
        DESCRIPTION    Return the count of code_name_secondary_structure pair for interaction type.
                       'code_name_secondary_structure' is set as index of the table to easly merge with other interaction type row (or table).
        RETURN         pandas dataframe
        """        
        # Count the occurrences of each element
        element_counts = Counter(self.list_interaction_code_simple)
        
        # Convert the Counter object to a Pandas DataFrame
        df = pd.DataFrame(list(element_counts.items()), columns=['pair_code_name_secondary_structure', self.result.__class__.__name__.lower()])
        df.set_index('pair_code_name_secondary_structure', inplace=True)
        
        return df

    
    
    @property
    def get_table_code_complete(self):
        """
        DESCRIPTION    Return the count of code_complete pair for interaction type.
                       'code_complete' is set as index of the table to easly merge with other interaction type row (or table).
        RETURN         pandas dataframe
        """        
        # Count the occurrences of each element
        element_counts = Counter(self.list_interaction_code_complete)
        
        # Convert the Counter object to a Pandas DataFrame
        df = pd.DataFrame(list(element_counts.items()), columns=['pair_code_complete', self.result.__class__.__name__.lower()])
        df.set_index('pair_code_complete', inplace=True)
        
        return df




#=====================================================
#===== END
#=====================================================
if __name__ == "__main__":
    main()
