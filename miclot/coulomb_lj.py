#!/usr/bin/env python3

"""
This script is part of MICLOT and is use to calculate Coulomb and Lennard-Jones energies
between residues.
"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"


__all__ = ['coulomb_lj','omm_coulomb_lj']


#=====================================================
#===== Import modules
#=====================================================

import xml.etree.ElementTree as ET
import mdtraj as md
import numpy as np

#openmm
from openmm import *
from openmm.app import *
from openmm.unit import *



#=====================================================
#===== Class for atom names conversion
#=====================================================

class convert_atom_name:
    def __init__(self, atom_name):
        """
        Convert atom names from AMBER to CHARMM
                        or from CHARMM to AMBER
        """
        #===== Initialise variable =====
        self.atom_name = atom_name
        
        #===== create convertion dictionnary AMBER to CHARMM =====
        self.conversion_dict = {
            'H'  : 'HN',
            'HB2': 'HB1',
            'HB3': 'HB2',
            'HD2': 'HD1',
            'HD3': 'HD1',
            'HE2': 'HE1',
            'HE3': 'HE2',
            'HG2': 'HG1',
            'HG3': 'HG2',
            'HZ2': 'HZ1',
            'HZ3': 'HZ2',
        }
        
        #===== create convertion dictionnary CHARMM to AMBER =====
        self.reverse_conversion_dict = {v: k for k, v in self.conversion_dict.items()}
        
        
        
    #===== Return parameters =====
    @property
    def amber_to_charmm(self):
        """
        DESCRIPTION    Convert AMBER atom name to CHARM atom name.
        """
        if self.atom_name in self.conversion_dict:
            return self.conversion_dict[self.atom_name]
        else:
            return self.atom_name

    @property
    def charmm_to_amber(self):
        """
        DESCRIPTION    Convert CHARM atom name to AMBER atom name.
        """        
        if self.atom_name in self.reverse_conversion_dict:
            return self.reverse_conversion_dict[self.atom_name]
        else:
            return self.atom_name





#=====================================================
#===== Class to parse CHARMM force field
#=====================================================

class charmm_atom_parameters:
    def __init__(self, forcefield_xml, residue_name, atom_name):
        """
        DEBUGGING
            AttributeError: 'NoneType' object has no attribute 'get'    It occure when the atom_name (or residue_name)
                                                                        don't exist in the Force Field. 
        """
        #===== Initialise variable =====
        self.forcefield = forcefield_xml
        self.residue_name = residue_name
        self.atom_name = atom_name
        #
        self.atom = self.forcefield.find(f".//Residue[@name='{self.residue_name}']/Atom[@name='{self.atom_name}']")
        self.atom_charge = self.atom.get('charge')
        self.atom_type = self.atom.get('type')
        self.atom_class = self.forcefield.find(f".//AtomTypes/Type[@name='{self.atom_type}']").get('class')

        #
        self.atom_LJForce = self.forcefield.find(f".//LennardJonesForce/Atom[@type='{self.atom_class}']")
        self.atom_LJForce_sigma = self.atom_LJForce.get('sigma')
        self.atom_LJForce_epsilon = self.atom_LJForce.get('epsilon')

    @property
    def charge(self):
        return self.atom_charge
    
    @property
    def sigma(self):
        return self.atom_LJForce_sigma
    
    @property
    def epsilon(self):
        return self.atom_LJForce_epsilon
    





#=====================================================
#===== Class to parse AMBER force field
#=====================================================

class amber_atom_parameters:
    def __init__(self, forcefield_xml, residue_name, atom_name):
        """ 
        DEBUGGING
            AttributeError: 'NoneType' object has no attribute 'get'    It occure when the atom_name (or residue_name)
                                                                        don't exist in the Force Field. 
        """
        #===== Initialise variable =====
        self.forcefield = forcefield_xml
        self.residue_name = residue_name
        self.atom_name = atom_name
        
        #
        self.atom = self.forcefield.find(f".//Residue[@name='{self.residue_name}']/Atom[@name='{self.atom_name}']")
        self.atom_charge = self.atom.get('charge')
        self.atom_type = self.atom.get('type')
        #self.atom_class = self.forcefield.find(f".//AtomTypes/Type[@name='{self.atom_type}']").get('class')

        #
        self.atom_LJForce = self.forcefield.find(f".//NonbondedForce/Atom[@type='{self.atom_type}']")
        self.atom_LJForce_sigma = self.atom_LJForce.get('sigma')
        self.atom_LJForce_epsilon = self.atom_LJForce.get('epsilon')

        
    #===== Return parameters =====
    @property
    def charge(self):
        return self.atom_charge
    
    @property
    def sigma(self):
        return self.atom_LJForce_sigma
    
    @property
    def epsilon(self):
        return self.atom_LJForce_epsilon





#=====================================================
#===== Class for Coulomb energy between 2 atoms
#=====================================================

class atoms_coulomb:
    def __init__(self, distance_A_B, charge_atom_A, charge_atom_B, solute_Dielectric, solvent_Dielectric):
        """
        DEBUGGING
            AttributeError: 'NoneType' object has no attribute 'get'    It occure when the atom_name (or residue_name)
                                                                        don't exist in the Force Field. 
        """
        #===== Initialise variable =====
        self.distance_A_B = float(distance_A_B)
        self.charge_atom_A = float(charge_atom_A)
        self.charge_atom_B = float(charge_atom_B)
        self.solute_Dielectric = float(solute_Dielectric)
        self.solvent_Dielectric = float(solvent_Dielectric)
        
        
        
    #===== Return parameters =====
    @property
    def get_energy(self):
        """
        DESCRIPTION    coulomb interaction without Cutoff.
                       All atom pair are considered.
        RETURN    energy
        UNIT      kJ/mol
        """
        # electric conversion factor f = 1/(4*Pi*epsilon_0) is set to 138.935458 kJ.mol-1.nm.e-2 
        return 138.935458 * ((self.charge_atom_A * self.charge_atom_B) / (self.distance_A_B * self.solute_Dielectric)) 

    
    #-----
    def cutoff(self,cutoff):
        """
        DESCRIPTION    coulomb interaction with 'hard' Cutoff.
                       All distance greater than the cutoff are ignored.
                       
        ARGUMENTS
            cutoff    Distance in nm.
            
        RETURN    energy
        UNIT      kJ/mol
        """
        if self.distance_A_B <= cutoff:
            # electric conversion factor f = 1/(4*Pi*epsilon_0) is set to 138.935458 kJ.mol-1.nm.e-2 
            return 138.935458 * ((self.charge_atom_A * self.charge_atom_B) / (self.distance_A_B * self.solute_Dielectric))
        else:
            return 0
   

    #-----
    def cutoff_reactionfield(self,cutoff):
        """
        DESCRIPTION    coulomb interaction with Cutoff.
                       Use reaction field.
        
        ARGUMENTS
            cutoff    Distance in nm.
            
        RETURN    energy
        UNIT      kJ/mol
        """
        k_rf = (1 / cutoff**3) * ((self.solvent_Dielectric - self.solute_Dielectric) / (2 * self.solvent_Dielectric + self.solute_Dielectric))
        c_rf = (1 / cutoff) * ((3*self.solvent_Dielectric) / (2*self.solvent_Dielectric + self.solute_Dielectric))
        
        # electric conversion factor f = 1/(4*Pi*epsilon_0) is set to 138.935458 kJ.mol-1.nm.e-2
        return 138.935458 * (self.charge_atom_A*self.charge_atom_B /self.solute_Dielectric) * ((1/self.distance_A_B) + (k_rf * self.distance_A_B**2) - c_rf)

        



#=====================================================
#===== Class for LJ energy between 2 atoms
#=====================================================

class atoms_lennard_jones:
    def __init__(self, distance_A_B, sigma_atom_A, sigma_atom_B, epsilon_atom_A, epsilon_atom_B):
        """
        DEBUGGING
            AttributeError: 'NoneType' object has no attribute 'get'    It occure when the atom_name (or residue_name)
                                                                        don't exist in the Force Field. 
        """
        #===== Initialise variable =====
        self.distance_A_B = float(distance_A_B)
        self.sigma_atom_A = float(sigma_atom_A)
        self.sigma_atom_B = float(sigma_atom_B)
        self.epsilon_atom_A = float(epsilon_atom_A)
        self.epsilon_atom_B = float(epsilon_atom_B)
    
        #===== Compute energy =====
        self.sigma = (self.sigma_atom_A + self.sigma_atom_B)/2
        self.epsilon = np.sqrt([self.epsilon_atom_A * self.epsilon_atom_B])[0]
        self.energy = 4 * self.epsilon * ( (self.sigma/self.distance_A_B)**12 - (self.sigma/self.distance_A_B)**6 )
    
    
    
    #===== Return parameters =====
    @property    
    def get_energy(self):
        """
        DESCRIPTION    Lennard-Jones interaction between each pair of atoms.
        RETURN    energy
        UNIT      kJ/mol
        """
        return self.energy
    
    #-----
    def switching_fucntion(self, d_switch, cutoff):
        """
        DESCRIPTION    Lennard-Jones interaction between each pair of atoms,
                       using a switching function to make the energy go
                       slowly to 0 at the cutoff distance.
                       
        ARGUMENT
            cutoff    Distance in nm. 
            
        RETURN    energy
        UNIT      kJ/mol
        """
        if cutoff <= d_switch:
            raise ValueError("Switching must be lower than the cutoff.")
            
        elif d_switch <= self.distance_A_B and self.distance_A_B < cutoff:
            x = (self.distance_A_B - d_switch) / (cutoff - d_switch)
            s = 1 - (6 * x**5) + (15 * x**4) - (10 * x**3)
            return s * self.energy
        
        elif self.distance_A_B < d_switch:
            return self.energy
        
        elif cutoff <= self.distance_A_B:
            return 0
    
    #-----
    def cutoff(self,cutoff):
        """
        DESCRIPTION    Lennard-Jones interaction with 'hard' Cutoff.
                       All distance greater than the cutoff are ignored.
                       
        ARGUMENTS
            cutoff    Distance in nm.
            
        RETURN    energy
        UNIT      kJ/mol
        """
        if self.distance_A_B <= cutoff: 
            return self.energy
        else:
            return 0
        
        
        
    @property    
    def get_sigma(self):
        """
        DESCRIPTION    Return the sigma of the given atom pair.
                       sigma = (sigma_A + sigma_B) / 2
        """
        return self.sigma

    
    @property    
    def get_epsilon(self):
        """
        DESCRIPTION    Return the epsilon of the given atom pair.
                       epsilon = sqrt( epsilon_A * epsilon_B )
        """
        return self.epsilon
    




#=====================================================
#===== Class to calculate Coulomb & LJ energies between 2 residues
#=====================================================

class coulomb_lj:
    def __init__(self, trajectory, res_index_A, res_index_B, force_field_path, frame=0, coulomb_cutoff=12, lj_cutoff=12, \
                 lj_switch=8, solute_dielectric=1.0, solvent_dielectric=78.5):
        """
        ENERGY TYPE    Coulomb & Lennar-Jones
        SUBTYPE(S)     Coulomb, Lennar-Jones

        DESCRIPTION
            Coulomb & Lennar-Jones interaction between residues.
            
        ARGUMENTS
            trajectory     MDTraj trajectory
            res_index_A    Index of residue A
            res_index_B    Index of residue B
        
        OPTIONAL ARGUMENTS
            frame                 Frame ID on which to perform the analysis.
                                  Default value: 0
                                  
            force_field_path      Path of the forcefiel file to use.
                                  For AMBER use hte file:  protein.ff14SB.xml
                                  For CHARMM use the file: charmm36.xml
                                  
            coulomb_cutoff        Cutoff applied to coulombic energy.
                                  It is used for hard cutoff and for reaction field.
                                  Unit: Å
                                  Default value: 12
                                  
            lj_cutoff             Cutoff applied to LJ energy.
                                  It is used for hard cutoff and for switching function.
                                  Unit: Å
                                  Default value: 12
                                  
            lj_switch             Switch distance applied to LJ energy.
                                  Unit: Å
                                  Default value: 8
                                  
            solute_dielectric     solute dielectric constant
                                  Default value: 1.0
                                  
            solvent_dielectric    solvent dielectric constant
                                  Default value: 78.5
        """
        #===== Initialise variables =====
        self.traj = trajectory[frame]
        self.top = self.traj.topology
        self.force_field_path = force_field_path
        self.solute_dielectric = solute_dielectric
        self.solvent_dielectric = solvent_dielectric
        # /10 to convert all distance values from Å to nm
        self.coulomb_cutoff = coulomb_cutoff/10
        self.lj_cutoff = lj_cutoff/10
        self.lj_switch = lj_switch/10
        # Residues
        self.res_index_A = res_index_A
        self.res_index_B = res_index_B
        self.residue_A = self.top.residue(self.res_index_A)
        self.residue_B = self.top.residue(self.res_index_B)
        
        
        #===== Initialize energy variables =====
        # Lennard-Jones
        self.energy_LJ = 0
        self.energy_LJ_switch = 0
        self.energy_LJ_cutoff = 0
        # Coulomb
        self.energy_C = 0
        self.energy_C_cutoff = 0
        self.energy_C_reactionfield = 0
        
        
        #===== Parse force field xml file =====
        self.ffxml = ET.parse(self.force_field_path)
        self.force_field = self.ffxml.getroot()
        
        
        #===== Calculate Coulomb and LJ energies =====
        # loop over all atoms in residue A
        for self.atom_A in self.residue_A.atoms:
            
            # get convert atom name to correspond to selected force field
            if 'charmm36.xml' in self.force_field_path:
                self.atom_A_name = convert_atom_name(self.atom_A.name).amber_to_charmm
                # atoms proterties from the force field
                self.atom_A_properties = charmm_atom_parameters(self.force_field, self.residue_A.name, self.atom_A_name)
            elif 'protein.ff14SB.xml' in self.force_field_path:
                self.atom_A_name = convert_atom_name(self.atom_A.name).charmm_to_amber
                self.atom_A_properties = amber_atom_parameters(self.force_field, self.residue_A.name, self.atom_A_name)
            else:
                raise ValueError('Force field file name must be charmm36.xml or protein.ff14SB.xml')

            
            # loop over all atoms in residue B
            for self.atom_B in self.residue_B.atoms:

                # get convert atom name to correspond to selected force field
                if 'charmm36.xml' in self.force_field_path:
                    self.atom_B_name = convert_atom_name(self.atom_B.name).amber_to_charmm
                    # atoms proterties from the force field
                    self.atom_B_properties = charmm_atom_parameters(self.force_field, self.residue_B.name, self.atom_B_name)
                elif 'protein.ff14SB.xml' in self.force_field_path:
                    self.atom_B_name = convert_atom_name(self.atom_B.name).charmm_to_amber
                    # atoms proterties from the force field
                    self.atom_B_properties = amber_atom_parameters(self.force_field, self.residue_B.name, self.atom_B_name)
                else:
                    raise ValueError('Force field file name must be charmm36.xml or protein.ff14SB.xml')

                
                # measure distance betwee atom_A and atom_B
                self.distance = md.compute_distances(self.traj, [[self.atom_A.index, self.atom_B.index]])[0][0]
                
                # calculate Coulomb and lennard-jones energie for the atom pair
                self.interaction_C = atoms_coulomb(self.distance, self.atom_A_properties.charge, self.atom_B_properties.charge, solute_Dielectric=self.solute_dielectric, solvent_Dielectric=self.solvent_dielectric)
                self.interaction_LJ = atoms_lennard_jones(self.distance, self.atom_A_properties.sigma, self.atom_B_properties.sigma, self.atom_A_properties.epsilon, self.atom_B_properties.epsilon)

                # Update LJ energy for the residue pair
                self.energy_LJ += self.interaction_LJ.get_energy
                self.energy_LJ_cutoff += self.interaction_LJ.cutoff(self.lj_cutoff)
                self.energy_LJ_switch += self.interaction_LJ.switching_fucntion(d_switch=self.lj_switch, cutoff=self.lj_cutoff)
                
                # Update Coulomb energy for the residue pair
                self.energy_C += self.interaction_C.get_energy
                self.energy_C_cutoff += self.interaction_C.cutoff(self.coulomb_cutoff)
                self.energy_C_reactionfield += self.interaction_C.cutoff_reactionfield(self.coulomb_cutoff)

                
        
    #===== Return results =====
    @property
    def get_energy(self):
        """
        DESCRIPTION    Return the total energy of the pair (LJ + Coulomb)
                       without considering cutoff and with the "hard" cutoff for LJ and Coulomb.
        RETURN         total, total_cutoff
        UNIT           kJ/mol
        """
        total = self.energy_LJ + self.energy_C
        total_cutoff = self.energy_LJ_cutoff + self.energy_C_cutoff
        return total, total_cutoff
    
    @property
    def get_energy_LJ(self):
        """
        DESCRIPTION    Return all calculated energy for LJ.
        RETURN         energy_LJ, energy_LJ_cutoff, energy_LJ_switch
        UNIT           kJ/mol
        """
        return self.energy_LJ, self.energy_LJ_cutoff, self.energy_LJ_switch
    
    @property
    def get_energy_coulomb(self):
        """
        DESCRIPTION    Return all calculated energy for Coulomb.
        RETURN         energy_C, energy_C_cutoff, energy_C_reactionfield
        UNIT           kJ/mol
        """
        return self.energy_C, self.energy_C_cutoff, self.energy_C_reactionfield    
    




#=====================================================
#===== Class to calculate Coulomb and LJ using OpenMM
#=====================================================
class omm_coulomb_lj:
    def __init__(self, trajectory, index_residue_A, index_residue_B, method=NoCutoff, nonbonded_cutoff=10.0, frame=0):
        """
        DESCRIPTION
            Adapted from https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/Computing%20Interaction%20Energies.html
            
            Calculate Coulomb and Lennard-Jones energies using openMM.
            It automaticaly return values for AMBER and CHARMM force fields.
        
        
        ARGUMENTS
            trajectory         MDTraj trajectory
            index_residue_A    MDTraj index of residue A
            index_residue_B    MDTraj index of residue B
            
            
        OPTIANL ARGUMENTS
            frame     MDTraj frame ID on which to perform the analysis.
            
            method    openMM keyword (not a string) defining the method to use for nonbonded interactions.
                      Allowed values: NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, or LJPME
                      Default value: NoCutoff
                      
            nonbonded_cutoff    cutoff distance to use for nonbonded interactions
                                Default value: 10 angstrom
        """        
        #===== Initialise variable =====
        self.trajectory = trajectory[frame]
        self.index_residue_A = index_residue_A
        self.index_residue_B = index_residue_B
        self.method = method
        self.nonbonded_cutoff = nonbonded_cutoff /10 # /10 to convert angstrom into nm
        
        #===== Initialize platfom and platform properties used by openMM =====
        self.platform = Platform.getPlatformByName('CPU')     # force using only CPU. Using CUDA or OpenCL will slowdown due memory transfer.
        self.platform.setPropertyDefaultValue('Threads', '1') # force using only one thread (more will slowdown the speed)

        
        #===== Convert MDTraj topology and position to openMM =====
        #convert the topology
        self.openmm_topology = self.trajectory.topology.to_openmm()
        
        # Get the positions from the first frame (for example)
        self.positions = self.trajectory.xyz[0]  # MDTraj stores positions in nanometers

        # Convert MDTraj positions to OpenMM positions
        self.openmm_positions = [Vec3(pos[0], pos[1], pos[2]) *nanometer for pos in self.positions]
        
        
        #===== calculate energies =====
        self.amber_total, self.amber_coulomb, self.amber_lj = self.calculate_amber()
        self.charmm_total, self.charmm_coulomb, self.charmm_lj = self.calculate_charmm()
        
        
        
    #===== Functions ======
    #----- scaling CHARMM Coulomb------
    def scaling_charmm_energy_coulomb(self, context, residue_A_scale, residue_B_scale):
        """
        Extract the desired par in the energy. 0: inactivated / 1: activated
        """
        context.setParameter("residue_A_scale", residue_A_scale)
        context.setParameter("residue_B_scale", residue_B_scale)
        return context.getState(getEnergy=True, groups={0}).getPotentialEnergy()
    
    
    #----- scaling AMBER Coulomb and LJ -----
    def scaling_amber_energy(self, context, residue_A_coulomb_scale, residue_A_lj_scale, residue_B_coulomb_scale, residue_B_lj_scale):
        """
        Extract the desired par in the energy. 0: inactivated / 1: activated
        """
        context.setParameter("residue_A_coulomb_scale", residue_A_coulomb_scale)
        context.setParameter("residue_A_lj_scale", residue_A_lj_scale)
        context.setParameter("residue_B_coulomb_scale", residue_B_coulomb_scale)
        context.setParameter("residue_B_lj_scale", residue_B_lj_scale)
        return context.getState(getEnergy=True, groups={0}).getPotentialEnergy()
    
    
    #----- Calculate energy with CHARMM -----
    def calculate_charmm(self):
        """
        Calculate Coulomb and LJ energies using CHARMM force field.
        """
        # Create the system using charmm FF and with/without cutoff
        forcefield = ForceField('charmm36.xml')
        system = forcefield.createSystem(self.openmm_topology, nonbondedMethod=self.method, nonbondedCutoff=self.nonbonded_cutoff*nanometer)

        # Define the selections for residues pair
        residue_A = set([i.index for i in self.openmm_topology.atoms() if i.residue.index == self.index_residue_A])
        residue_B = set([i.index for i in self.openmm_topology.atoms() if i.residue.index == self.index_residue_B])
        
        # Modify forces
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                force.setForceGroup(0)
                force.addGlobalParameter("residue_A_scale", 1)
                force.addGlobalParameter("residue_B_scale", 1)
                
                for i in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(i)
                    # Set the parameters to be 0 when the corresponding parameter is 0,
                    # and to have their normal values when it is 1.
                    force.setParticleParameters(i, 0, 0, 0)
                    if i in residue_A:
                        force.addParticleParameterOffset("residue_A_scale", i, charge, sigma, epsilon)
                    elif i in residue_B:
                        force.addParticleParameterOffset("residue_B_scale", i, charge, sigma, epsilon)
                
                for i in range(force.getNumExceptions()):
                    p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
                    force.setExceptionParameters(i, p1, p2, 0, 0, 0)
            
            elif isinstance(force, CustomNonbondedForce):
                force.setForceGroup(1)
                force.addInteractionGroup(residue_A, residue_B)
            
            else:
                force.setForceGroup(2)


        # Create a Context for performing calculations.
        #     The integrator is not important, since we will only be performing single point energy evaluations.
        integrator = VerletIntegrator(0.001*picosecond)
        context = Context(system, integrator, self.platform)
        context.setPositions(self.openmm_positions)
        
        # Get the total coulomb energies for each residues
        total_coulomb     = self.scaling_charmm_energy_coulomb(context, 1, 1)
        residue_A_coulomb = self.scaling_charmm_energy_coulomb(context, 1, 0)
        residue_B_coulomb = self.scaling_charmm_energy_coulomb(context, 0, 1)
        
        # Calculate coulomb energy specificly for the pair
        coulomb_pair = total_coulomb - residue_A_coulomb - residue_B_coulomb
        
        # Calculate LJ energy specificly for the pair
        lj_pair = context.getState(getEnergy=True, groups={1}).getPotentialEnergy()
        
        # calculate coulomb + LJ for the pair
        lj_coulomb_pair = coulomb_pair + lj_pair

        # return the energies
        return lj_coulomb_pair, coulomb_pair, lj_pair
    
    
    #----- Calculate energy with AMBER -----------
    def calculate_amber(self):
        """
        Calculate Coulomb and LJ energies using AMBER force field.
        """
        # Create the system using amber FF and with/without cutoff
        forcefield = ForceField('amber14-all.xml')
        system = forcefield.createSystem(self.openmm_topology, nonbondedMethod=self.method, nonbondedCutoff=self.nonbonded_cutoff*nanometer)

        # Define the selections for residues pair
        residue_A = set([i.index for i in self.openmm_topology.atoms() if i.residue.index == self.index_residue_A])
        residue_B = set([i.index for i in self.openmm_topology.atoms() if i.residue.index == self.index_residue_B])
        
        # Modify forces
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                force.setForceGroup(0)
                force.addGlobalParameter("residue_A_coulomb_scale", 1)
                force.addGlobalParameter("residue_A_lj_scale", 1)
                force.addGlobalParameter("residue_B_coulomb_scale", 1)
                force.addGlobalParameter("residue_B_lj_scale", 1)
                for i in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(i)
                    force.setParticleParameters(i, 0, 0, 0)
                    if i in residue_A:
                        force.addParticleParameterOffset("residue_A_coulomb_scale", i, charge, 0, 0)
                        force.addParticleParameterOffset("residue_A_lj_scale", i, 0, sigma, epsilon)
                    elif i in residue_B:
                        force.addParticleParameterOffset("residue_B_coulomb_scale", i, charge, 0, 0)
                        force.addParticleParameterOffset("residue_B_lj_scale", i, 0, sigma, epsilon) 
                for i in range(force.getNumExceptions()):
                    p1, p2, chargeProd, sigma, epsilon = force.getExceptionParameters(i)
                    force.setExceptionParameters(i, p1, p2, 0, 0, 0)
            else:
                force.setForceGroup(2)

        ## Create a Context for performing calculations.
        ##     The integrator is not important, since we will only be performing single point energy evaluations.
        integrator = VerletIntegrator(0.001*picosecond)
        context = Context(system, integrator, self.platform)
        context.setPositions(self.openmm_positions)
        
        # Get coulomb energy
        coulomb_residue_A = self.scaling_amber_energy(context, 1, 0, 0, 0)
        coulomb_residue_B = self.scaling_amber_energy(context, 0, 0, 1, 0)
        coulomb_total     = self.scaling_amber_energy(context, 1, 0, 1, 0)
        coulomb_pair      = coulomb_total - coulomb_residue_A - coulomb_residue_B
        
        # Get LJ energy
        lj_residue_A = self.scaling_amber_energy(context, 0, 1, 0, 0)
        lj_residue_B = self.scaling_amber_energy(context, 0, 0, 0, 1)
        lj_total     = self.scaling_amber_energy(context, 0, 1, 0, 1)
        lj_pair = lj_total - lj_residue_A - lj_residue_B
        
        # Get total (coulomb + LJ) energy
        lj_coulomb_residue_A = self.scaling_amber_energy(context, 1, 1, 0, 0)
        lj_coulomb_residue_B = self.scaling_amber_energy(context, 0, 0, 1, 1)
        lj_coulomb_total     = self.scaling_amber_energy(context, 1, 1, 1, 1)
        lj_coulomb_pair      = lj_coulomb_total - lj_coulomb_residue_A - lj_coulomb_residue_B

        ## return the resulting energies
        return lj_coulomb_pair, coulomb_pair, lj_pair
    
    
    
    #===== Return energies =====
    @property
    def get_energy(self):
        """
        DESCRIPTION    Return total energy as calculated using AMBER and CHARMM force field
        RETURN         amber_total, charmm_total
        UNIT           Kj/mol
        """
        return self.amber_total._value, self.charmm_total._value
    
    @property
    def get_energy_coulomb(self):
        """
        DESCRIPTION    Return Coulomb energy as calculated using AMBER and CHARMM force field
        RETURN         amber_coulomb, charmm_coulomb
        UNIT           Kj/mol
        """
        return self.amber_coulomb._value, self.charmm_coulomb._value
    
    @property
    def get_energy_LJ(self):
        """
        DESCRIPTION    Return Coulomb energy as calculated using AMBER and CHARMM force field
        RETURN         amber_lj, charmm_lj
        UNIT           Kj/mol
        """
        return self.amber_lj._value, self.charmm_lj._value
        
        





#=====================================================
#===== END
#=====================================================
if __name__ == "__main__":
    main()
