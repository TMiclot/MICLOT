#!/usr/bin/env python3

"""
This script is part of MICLOT ...
"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"





#=====================================================
#===== Import modules
#=====================================================

import numpy as np
import mdtraj as md
import pandas as pd
from tqdm import tqdm





#=====================================================
#===== Class for bridge geometric parameters
#=====================================================

class bridge_geomerty:
    def __init__(self, trajectory, residue_index_A, residue_index_B):
        """
        DESCRIPTION
            Return geometric parameters between two Cys:
                - X1, X2, X3 and theta angles
                - CA, CB, XX distances (where X is element S or Se)
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of cysteine A
            res_index_B    Index of cysteine B
        """
        #===== Initialise variables =====
        self.traj = trajectory
        self.top = self.traj.topology
        self.residue_A = residue_index_A
        self.residue_B = residue_index_B
        
        #===== Get atoms index for residue A ===== 
        self.atom_CA_resA = self.top.select(f"resid {self.residue_A} and name CA")[0]
        self.atom_CB_resA = self.top.select(f"resid {self.residue_A} and name CB")[0]
        self.atom_N_resA = self.top.select(f"resid {self.residue_A} and name N")[0]
        self.atom_XX_resA = self.top.select(f"resid {self.residue_A} and (name SG or name SE)")[0] # XX is sulfur or selenium
        
        #===== Get atoms index for residue A ===== 
        self.atom_CA_resB = self.top.select(f"resid {self.residue_B} and name CA")[0]
        self.atom_CB_resB = self.top.select(f"resid {self.residue_B} and name CB")[0]
        self.atom_N_resB = self.top.select(f"resid {self.residue_B} and name N")[0]
        self.atom_XX_resB = self.top.select(f"resid {self.residue_B} and (name SG or name SE)")[0] # XX is sulfur or selenium

        
    #==== DIHEDRALS ====
    @property
    def angle_X1(self):
        """
        DIHEDRAL ANGLE    N-CA-CB-XX (XX is S or Se)
        RETURN            Angle value.
                          Values for residue 1 and residue 2.
        UNIT              degree
        """
        # set atom indices
        atoms_indices_A = [[self.atom_N_resA, self.atom_CA_resA, self.atom_CB_resA, self.atom_XX_resA]]
        atoms_indices_B = [[self.atom_N_resB, self.atom_CA_resB, self.atom_CB_resB, self.atom_XX_resB]]
        # calculate angles
        angle_A = md.compute_dihedrals(self.traj, atoms_indices_A)[0][0]
        angle_B = md.compute_dihedrals(self.traj, atoms_indices_B)[0][0]
        # convert angle from radian to degree
        angle_A = np.rad2deg(angle_A)
        angle_B = np.rad2deg(angle_B)
        return angle_A, angle_B
        
    @property
    def angle_X2(self):
        """
        DIHEDRAL ANGLE    CA-CB-XX1-XX2 (XX is S or Se)
        RETURN            Angle value.
                          Values for residue 1 and residue 2.
        UNIT              degree
        """
        # set atom indices
        atoms_indices_A = [[self.atom_CA_resA, self.atom_CB_resA, self.atom_XX_resA, self.atom_XX_resB]]
        atoms_indices_B = [[self.atom_CA_resB, self.atom_CB_resB, self.atom_XX_resB, self.atom_XX_resA]]
        # calculate angles
        angle_A = md.compute_dihedrals(self.traj, atoms_indices_A)[0][0]
        angle_B = md.compute_dihedrals(self.traj, atoms_indices_B)[0][0]
        # convert angle from radian to degree
        angle_A = np.rad2deg(angle_A)
        angle_B = np.rad2deg(angle_B)
        return angle_A, angle_B
    
    @property
    def angle_X3(self):
        """
        DIHEDRAL ANGLE    CB2-XX1-XX2-CB2 (XX is S or Se)
        RETURN            Angle value.
        UNIT              degree
        """
        # set atom indices
        atoms_indices = [[self.atom_CB_resA, self.atom_XX_resA, self.atom_XX_resB, self.atom_CB_resB]]
        # measure the angle
        angle = md.compute_dihedrals(self.traj, atoms_indices)[0][0]
        # convert angle from radian to degree
        angle = np.degrees(angle)
        return angle
    
    
    #==== ANGLES ====
    @property
    def angle_theta(self):
        """
        ANGLE     CA-CB-XX (XX is S or Se)
        RETURN    Angle value.
                  Values for residue 1 and residue 2.
        UNIT      degree
        """
        # Set atom indices
        atoms_indices_A = [[self.atom_CA_resA, self.atom_CB_resA, self.atom_XX_resA]]
        atoms_indices_B = [[self.atom_CA_resB, self.atom_CB_resB, self.atom_XX_resB]]
        
        # measure the angles
        angle_A = md.compute_angles(self.traj, atoms_indices_A)[0][0]
        angle_B = md.compute_angles(self.traj, atoms_indices_B)[0][0]
        
        # convert angle from radian to degree
        angle_A = np.rad2deg(angle_A)
        angle_B = np.rad2deg(angle_B)
        return angle_A, angle_B
    
    
    #==== DISTANCES ====
    @property
    def distance_CA(self):
        """
        DISTANCE    CA-CA between the two residues.
        UNIT        Angstrom (A)
        """
        atom_pairs = [[self.atom_CA_resA, self.atom_CA_resB]]
        # *10 is used to convert nm to angstrom
        distance = md.compute_distances(self.traj, atom_pairs)[0][0] *10
        return distance
    
    @property
    def distance_CB(self):
        """
        DISTANCE    CB-CB between the two residues.
        UNIT        Angstrom (A)
        """
        atom_pairs = [[self.atom_CB_resA, self.atom_CB_resB]]
        # *10 is used to convert nm to angstrom
        distance = md.compute_distances(self.traj, atom_pairs)[0][0] *10
        return distance
    
    @property
    def distance_X(self):
        """
        DISTANCE    S-S or Se-Se or S-Se between the two residues.
        UNIT        Angstrom (A).
        """
        atom_pairs = [[self.atom_XX_resA, self.atom_XX_resB]]
        # *10 is used to convert nm to angstrom
        distance = md.compute_distances(self.traj, atom_pairs)[0][0] *10
        return distance





#=====================================================
#===== Class for bridge energies
#=====================================================

class bridge_energy_SS:
    def __init__(self, x1_A, x1_B, x2_A, x2_B, x3, ta, tb):
        """
        DESCRIPTION
            Return the energy of a Cys-Cys bridge, only for disulfide type.
            It can calculate the energy using 2 methods:
                - Disulfide model structure (equivalent to Disulfid by Design)
                - Dihedral strain energy 
        
        ARGUMENTS
            x1_A    X1 angle of cysteine A
            x1_B    X1 angle of cysteine B
            x2_A    X2 angle of cysteine A
            x2_B    X2 angle of cysteine B
            x3      X3 angle between cysteines A and B
            ta      theta angle of cysteine A
            tb      theta angle of cysteine B
            
            * angles are in degree        
        """
        #===== Initialise variables =====
        self.x1_A = x1_A # Angle X1 of residue A
        self.x1_B = x1_B # Angle X1 of residue B
        self.x2_A = x2_A # Angle X2 of residue A
        self.x2_B = x2_B # Angle X2 of residue B
        self.x3 = x3   # Angle X3 between residues A and B
        self.t_a = ta   # Angle theta of residue A
        self.t_b = tb   # Angle theta of residue B
        
        
        #===== Compute: Disulfide model structure (DMS) Part =====
        self.E_dms_x1_A = 1.4*(1 + np.cos(np.deg2rad(3*self.x1_A)))
        self.E_dms_x1_B = 1.4*(1 + np.cos(np.deg2rad(3*self.x1_B)))
        
        self.E_dms_x3 = 4.0*(1 - np.cos(np.deg2rad(1.957*(self.x3 + 87.0))))
        
        self.E_dms_t_a = 55.0 * (self.t_a - 114.6)**2 /4184
        self.E_dms_t_b = 55.0 * (self.t_b - 114.6)**2 /4184
        
        
        #===== Compute: Dihedral strain energy (DSE) Part ====
        self.E_dse_x1_A = 8.37*(1 + np.cos(3*np.deg2rad(self.x1_A)))
        self.E_dse_x1_B = 8.37*(1 + np.cos(3*np.deg2rad(self.x1_B)))
        
        self.E_dse_x2_A = 4.18*(1 + np.cos(3*np.deg2rad(self.x2_A)))
        self.E_dse_x2_B = 4.18*(1 + np.cos(3*np.deg2rad(self.x2_B)))
        
        self.E_dse_x3 = 14.64*(1 + np.cos(2*np.deg2rad(self.x3))) + 2.51*(1 + np.cos(3*np.deg2rad(self.x3)))
    
    
    #===== Return dms results =====
    @property
    def get_dms_X1(self):
        """
        METHOD    Disulfide by Design v2
        ANGLE     N-CA-CB-S
        RETURN    Energy of the dihedral angle.
                  Values for residue 1 and residue 2.
        UNIT      kcal/mol
        """
        return self.E_dms_x1_A, self.E_dms_x1_B

    @property
    def get_dms_X3(self):
        """
        METHOD    Disulfide by Design v2
        ANGLE     CB1-XX1-XX2-CB2 (XX is S or Se)
        RETURN    Energy of the dihedral angle.
        UNIT      kcal/mol
        """
        return self.E_dms_x3

    @property
    def get_dms_theta(self):
        """
        METHOD    Disulfide by Design v2
        ANGLE     CA-CB-S
        RETURN    Energy of the angle.
                  Values for residue 1 and residue 2.
        UNIT      kcal/mol
        """
        return self.E_dms_t_a, self.E_dms_t_b

    @property
    def get_dms_total(self):
        """
        METHOD    Disulfide by Design v2
        RETURN    Energy of the disulfide bridge.
        UNIT      kcal/mol
        """
        return self.E_dms_x1_A + self.E_dms_x1_B + self.E_dms_x3 + self.E_dms_t_a + self.E_dms_t_b
    
    
    #===== Return DSE results =====
    @property
    def get_dse_X1(self):
        """
        METHOD    Dihedral strain energy
        ANGLE     N-CA-CB-S
        RETURN    Energy of the angle.
                  Values for residue 1 and residue 2.
        UNIT      kJ/mol
        """
        return self.E_dse_x1_A, self.E_dse_x1_B

    @property
    def get_dse_X2(self):
        """
        METHOD    Dihedral strain energy
        ANGLE     CA-CB-S1-S2.
        RETURN    Energy of the angle.
                  Values for residue 1 and residue 2.
        UNIT      kJ/mol
        """
        return self.E_dse_x2_A, self.E_dse_x2_B    
    
    @property
    def get_dse_X3(self):
        """
        METHOD    Dihedral strain energy
        ANGLE     CB1-XX1-XX2-CB2 (XX is S or Se)
        RETURN    Energy of the angle.
                  Values for residue 1 and residue 2.
        UNIT      kJ/mol
        """
        return self.E_dse_x3
        
    @property
    def get_dse_total(self):
        """
        METHOD    Dihedral strain energy
        RETURN    Energy of the disulfide bridge.
        UNIT      kJ/mol
        """
        return self.E_dse_x1_A + self.E_dse_x1_B + self.E_dse_x2_A + self.E_dse_x2_B + self.E_dse_x3





#=====================================================
#===== Class to identify Cys bridge
#=====================================================

class cys_bridge:
    def __init__(self, trajectory, residue_index_A, residue_index_B, frame=0, MAX_distance_XX=None, MAX_distance_CA=7.5, \
                 MIN_distance_CA=3.0, MAX_distance_CB=5.5, method=None, MAX_energy=75.0, MIN_energy=0.0):
        """
        DESCRIPTION
            Identify disulfide, diselenide or selenosulfide bridge between two Cys.
        
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS
            frame              Frame in the traj
                               Default value: 0
            
            MAX_distance_XX    Maximum distance beween S, Se, atoms.
                               Default value: covalent radii
            MAX_distance_CA    Maximum distance between the CA of each Cys.
                               Default value: 7.5
            MIN_distance_CA    Minimum distance between the CA of each Cys.
                               Default value: 3.0
            MAX_distance_CB    Maximum distance between the CB of each Cys.
                               Default value: 5.5
                               
            method             Method to calculate energy of the bridge.
                               Use this argument only if you want to add energetic
                               restrain to identify Cys bridge.
                               Values are: None, 'DMS', 'DSE'
                               Default value: None
            MAX_energy         Maximum energy to consider a bridge.
                               To use it, remember to choose a method.
                               Default value: 75.0
            MIN_energy         Maximum energy to consider a bridge.
                               To use it, remember to choose a method.
                               Default value: 0.0
        """
        #===== Initialise variable =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.residue_A = residue_index_A
        self.residue_B = residue_index_B
        self.method = method
        self.MAX_energy = MAX_energy
        self.MIN_energy = MIN_energy
        self.MAX_distance_CA = MAX_distance_CA
        self.MIN_distance_CA = MIN_distance_CA
        self.MAX_distance_CB = MAX_distance_CB
        self.MAX_distance_XX = MAX_distance_XX
        
        
        
        
        #===== Check type of bridge: Cys-S--S-Cys, Cys-Se--Se-Cys, Cys-Se--S-Cys =====
        self.resiude_A_name = self.traj.topology.residue(self.residue_A).name
        self.resiude_B_name = self.traj.topology.residue(self.residue_B).name
        
        if self.resiude_A_name == "CYS" and self.resiude_B_name == "CYS":
            self.MAX_distance_XX_default = 2.10
            self.bridge_type = "disulfide"
        
        elif self.resiude_A_name == "SEC" and self.resiude_B_name == "SEC":
            self.MAX_distance_XX_default = 2.41
            self.bridge_type = "diselenide"
            
        elif self.resiude_A_name != self.resiude_B_name and self.resiude_A_name in ["CYS","SEC"] and self.resiude_B_name in ["CYS","SEC"]:
            self.MAX_distance_XX_default = 2.27
            self.bridge_type = "selenosulfide"
            
        
        
        #===== Check if user provide custom 'MAX_distance_XX', else take de default one =====
        if self.MAX_distance_XX == None:
            self.MAX_distance_XX = self.MAX_distance_XX_default
        
        
        #===== Check the interaction using distance parameters =====
        # compute geometric parameters
        self.geometry = bridge_geomerty(traj, self.residue_A, self.residue_B)

        # Geometric (distance) checking
        if (self.MIN_distance_CA <= self.geometry.distance_CA) and (self.geometry.distance_CA <= self.MAX_distance_CA) \
        and (self.geometry.distance_CB <= self.MAX_distance_CB) and (self.geometry.distance_X <= self.MAX_distance_XX):
            
            # Successfull pass the geometric (distance) checking
            self.geometric_checking = True

            # get angular parameters
            self.x1_i, self.x1_j = self.geometry.angle_X1
            self.x2_i, self.x2_j = self.geometry.angle_X2
            self.x3_ij = self.geometry.angle_X3
            self.theta_i, self.theta_j = self.geometry.angle_theta

            # compute energetic parameters ONLY for disulfide
            if self.bridge_type == "disulfide":
                try:
                    self.energies = bridge_energy_SS(self.x1_i, self.x1_j, self.x2_i, self.x2_j, self.x3_ij, self.theta_i, self.theta_j)
                except:
                    self.energies = None
            # if the bridge type is dielenide or selenosulfide, energy is not calcualted
            else:
                self.energies = None
        
        # Don't pass the geometric (distance) checking
        else:
            self.geometric_checking = False
            self.energies = None
 
           
        
    #===== Return results =====
    @property
    def check_interaction(self):
        """
        DESCRIPTION
            Check if two cyctein are invloved in a bridge (disulfide, diselenide, selenosulfide).
        
        RETURN
            True      The interaction exist.
            False     The interaction don't exist.
        """
        if self.method == None and self.geometric_checking == True:
            return True, self.bridge_type
        
        elif self.method == "DSE" and self.geometric_checking == True and (self.MIN_energy <= self.energies.get_dse_total) and (self.energies.get_dse_total <= self.MAX_energy):
            return True, self.bridge_type
        
        elif self.method == "DMS" and self.geometric_checking == True and (self.MIN_energy <= self.energies.get_dms_total) and (self.energies.get_dms_total <= self.MAX_energy):
            return True, self.bridge_type
        
        else:
            return False, False


    @property
    def get_energy(self):  
        """
        DESCRIPTION    Return the calculated energy for the method DSE and for the method DMS.
                       Warning: Their units are diffrent.
        UNIT           DSE: kJ/mol, DMS: kcal/mol
        RETURN         dse_total, dms_total
        """
        if self.energies != None:
            return self.energies.get_dse_total, self.energies.get_dms_total
        else:
            return None
    
    
    @property
    def get_energy_dse(self):
        """
        DESCRIPTION    Return the calculated energies for: X1 angles, X2 angles and X3 angle
                       used in the method DSE.
        UNIT           kJ/mol
        RETURN         (X1_A,  X1_B), (X2_A,  X2_B), X3
        """
        if self.energies != None:
            return self.energies.get_dse_X1, self.energies.get_dse_X2, self.energies.get_dse_X3
        else:
            return None
        
        
    @property
    def get_energy_dms(self):
        """
        DESCRIPTION    Return the calculated energies for: X1 angles, X3 angle and theta angles
                       used in the method DSE.
        UNIT           kcal/mol
        RETURN         (X1_A,  X1_B), X3, (theta_A,  theta_B)
        """
        if self.energies != None:
            return self.energies.get_dms_X1, self.energies.get_dms_X3, self.energies.get_dms_theta
        else:
            return None
        
        
    @property
    def get_distance(self):            
        """
        DESCRIPTION    Return all CA-CA, CB-CB and X-X distances.
        UNIT           angstrom
        RETURN         (X1_A,  X1_B), (X2_A,  X2_B), X3, (theta_a, theta_B)
        """
        return self.geometry.distance_CA, self.geometry.distance_CB, self.geometry.distance_X
        
        
    @property
    def get_angle(self):
        """
        DESCRIPTION    Return all X1, X2, X3 and theta angles.
        UNIT           Degree
        RETURN         (X1_A,  X1_B), (X2_A,  X2_B), X3, (theta_a, theta_B)
        """
        return self.geometry.angle_X1, self.geometry.angle_X2, self.geometry.angle_X3, self.geometry.angle_theta






#=====================================================
#===== Function to check bridge into an entire PDB file
#=====================================================

def check_cys_bridges(pdb_file, outfile=True, logfile=True):
    """
    DESCRIPTION
        This command check a PDB structure to identify all CYS-CYS bridges.
        It can identify disulfide, diselenide or selenosulfide bridges. 
        
    ARGUMENTS
        pdb_file    Protein structure in PDB format.
        
    OPTIONAL ARGUMENTS
        outfile     Generate a modifyes structure named '{NAME}_CYS_bridge.pdb'
                    with modified CYS to CYX, or SEC to XSE, names.
                    If set to False, warns the user that no PDB file is being written. 
                    Default value: True
        logfile     Generate an output file in CSV format containing all tests
                    performed and their results.
                    If set to False the result will be print in the terminal,
                    but with less informations.
                    Default value: True
    """
    #===== load PDB file in mdtraj =====
    traj = md.load(pdb_file, top=pdb_file)

    #===== Get index of all CYS and SEC in the structure (without redundancy) =====
    all_cys = set(traj.topology.atom(i).residue.index for i in traj.topology.select("protein and element S Se"))

    #===== Create a list of cys involved in bridge, based on distance only. =====
    list_bridges = []
    
    #===== Calculate total number of iterations =====
    total_iterations = len(all_cys) ** 2

    #===== Initialize tqdm with total number of iterations =====
    progress_bar = tqdm(total=total_iterations, desc='Cys-Cys Bridge - Progress')
    
    #===== If the logfile is True, create a pandas table to store all output file =====
    if logfile == True:
        df = pd.DataFrame(columns=['ResName_A', 'ResID_A', 'ResSeq_A', 'ChainID_A', 'ChainName_A',
                                   'ResName_B', 'ResID_B', 'ResSeq_B', 'ChainID_B', 'ChainName_B',
                                   'Bridge','Type'
                                  ])
    else:
        print("INFO - No logfile will be written.")

        
    #===== loop over all possible (seleno)cystein pairs =====
    for cys_A in all_cys:
        for cys_B in all_cys:
            
            #----- Update progress bar -----
            progress_bar.update(1)

            #----- Check the Cys-Cys bridge if the two Cys are not already involved into Cys Bridge -----
            if cys_A != cys_B and cys_A not in list_bridges and cys_B not in list_bridges:
                
                #----- Perform the checking -----
                bridge = cys_bridge(traj, cys_A, cys_B)

                #----- Rename Cys and get bridge type -----
                # If the bridge is disulfide, rename the two CYS to CYX
                if bridge.check_interaction[0] == True and bridge.check_interaction[1] == "disulfide":
                    # Add cys_A and cys_B to list_bridges
                    list_bridges.append(cys_A)
                    list_bridges.append(cys_B)
                    # Rename the two Cys in the topology
                    traj.topology.residue(cys_A).name = "CYX"
                    traj.topology.residue(cys_B).name = "CYX"

                # If the bridge is diselenide, rename the two CYS to XSE
                elif bridge.check_interaction[0] == True and bridge.check_interaction[1] == "diselenide":
                    # Add cys_A and cys_B to list_bridges
                    list_bridges.append(cys_A)
                    list_bridges.append(cys_B)
                    # Rename the two seleno-Cys in the topology
                    traj.topology.residue(cys_A).name = "XSE"
                    traj.topology.residue(cys_B).name = "XSE"

                # If the bridge is disulfide, rename the two CYS to CYX and SEC to XSE
                elif bridge.check_interaction[0] == True and bridge.check_interaction[1] == "selenosulfide":
                    # Add cys_A and cys_B to list_bridges
                    list_bridges.append(cys_A)
                    list_bridges.append(cys_B)
                    # Identify the Cys index and the seleno-Cys index
                    resid_XSE = traj.topology.select(f'resid {cys_A} {cys_B} and element Se')
                    resid_CYX = traj.topology.select(f'resid {cys_A} {cys_B} and element S')
                    # Rename the two (seleno-)Cys in the topology
                    traj.topology.residue(resid_XSE[0]).name = "XSE"
                    traj.topology.residue(resid_CYX[0]).name = "CYX"
                
                
                # if logfile is set to true, expand the pandas dataframe
                if logfile == True:
                    df_temp = pd.DataFrame({'ResName_A': [traj.topology.residue(cys_A).name],
                        'ResID_A': [cys_A],
                        'ResSeq_A': [traj.topology.residue(cys_A).resSeq],
                        'ChainID_A': [traj.topology.residue(cys_A).chain.index],
                        'ChainName_A': [chr(65+traj.topology.residue(cys_A).chain.index)],
                        'ResName_B': [traj.topology.residue(cys_B).name],
                        'ResID_B': [cys_B],
                        'ResSeq_B': [traj.topology.residue(cys_B).resSeq],
                        'ChainID_B': [traj.topology.residue(cys_B).chain.index],
                        'ChainName_B': [chr(65+traj.topology.residue(cys_B).chain.index)],
                        'Bridge': [str(bridge.check_interaction[0])],
                        'Type': [str(bridge.check_interaction[1])]
                        })
                    
                    df = pd.concat([df, df_temp], ignore_index=True)
                
                # if logfile is or set to true, print minimal information
                else:
                    print("\t[+]", cys_A, cys_B, bridge.check_interaction)
                    

    #===== Close the progress bar =====
    progress_bar.close()
    
    
    #===== Write the logfile =====
    if logfile == True:
        # Get the name of the pdb file
        name = pdb_file.split('.')[0]
        
        # Export DataFrame to CSV file
        df.to_csv(f'{name}_CYS_log.csv', index=False)  # Set index=False to exclude row indices

  
    #===== Export the the modified structure in PDB file or export the topology =====
    if outfile == True:
        # Get the name of the pdb file
        name = pdb_file.split('.')[0] 
        
        # Write PDB file frommodified topology
        traj.save(f"{name}_CYS_bridges.pdb")
        
    else:
        print("INFO - No structure files (.pdb) have been modified or written.")





#=====================================================
#===== END
#=====================================================
if __name__ == "__main__":
    main()
