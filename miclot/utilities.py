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
 
 This module is provide tools to prepare PDB files & to easily
 passthrough a PDB file and it MDTraj topology.

"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"





#=====================================================
#===== Import modules
#=====================================================

import os
import subprocess
#
import numpy as np
import mdtraj as md
import pandas as pd
# openMM
from openmm.app import *
from openmm import *
from openmm.unit import *





#=====================================================
#===== Function to convert chainID of MDTraj topology to chainName of PDB file
#=====================================================
def mdtraj_chainID_2_chainName(pdb_file_path, write_outfile=True):
    """
    DESCRIPTION
        Correlates a PDB file with its topology in MDTraj.

        The command returns two dictionaries:
        - convert chainID to chainName
        - convert chainName to chainID

        The command also writes three output files (as CSV):
        - correspondence between mdtraj toppology and PDB file
        - chainID and their corresponding chainName
        - chainName and their corresponding chainID

        Please note that depending on the PDB file, MDTraj may use different
        IDs for the same chain name in the PDB file.
        
        
    USAGE
        chainID_2_chainName, chainName_2_chainID = mdtraj_chainID_2_chainName('5azz.pdb')
    
    
    ARGUMENT
        pdb_file_path    Path to the pdb file (string format)
        
        
    OPTIONAL ARGUMENT
        write_outfile    (True/False) Write the output files.   
    """
    #===== Initialise column delimitation and names in PDB file for Pandas =====
    #----- PDB column delimitation -----
    # Columns: 12(11), 21(20), 28-30(27-29), 67-72(66-71) are empty in PDB format
    # SegmentID 73-75(72, 75) is not taked in account
    PDB_format_columns = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
                         (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (72, 76),
                         (76, 78), (78, 80)]
    
    #----- PDB column names -----
    PDB_format_contents = ['record', 'serial', 'name', 'alt_location', 'resName', 'chainName', 'resSeq',
                           'insertion', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'B_factor',
                           'segmentID', 'element', 'charge']
    
    
    
    #===== Read the PDB file with Pandas and get a Dataframe =====
    # Read a table of fixed-width formatted lines into DataFrame.
    df_pdb_file = pd.read_fwf(pdb_file_path, names=PDB_format_contents, colspecs=PDB_format_columns)

    # Get only ATOM and HETATM & only the 'alternative location' of 'A' or empty
    df_pdb_file_clean = df_pdb_file[(df_pdb_file['record'].isin(['ATOM', 'HETATM'])) & (df_pdb_file['alt_location'].isin(['A', np.NaN]))]
    
    # Reindex the dataframe
    df_pdb_file_clean.reset_index(drop=True, inplace=True)

    # Convers all values in the dataframe to string
    df_pdb_file_clean_string = df_pdb_file_clean.astype(str)
    
    
    
    #===== Read the PDB file with MDTraj and get a Dataframe =====
    # Get the PDB file as a topology
    topology = md.load(pdb_file_path, top=pdb_file_path).topology

    # Convert topology to dataframe
    df_mdtraj_residues, df_mdtraj_bonds = topology.to_dataframe()

    # Convert all values in the mdtraj_table to string
    df_mdtraj_residues_string = df_mdtraj_residues.astype(str)

    
    
    #===== Concatenate the MDTraj and PDB dataframes into one and save it as CSV =====
    # create a new dataframe with column of the 'df_pdb_file_clean_string' to keep and add to 'df_mdtraj_residues_string'
    extracted_columns = df_pdb_file_clean_string[['record', 'chainName', 'insertion', 'x_coord', 'y_coord', 'z_coord', 'occupancy', 'B_factor', 'charge']]
    
    # concatenate the mdtraj topology dataframe with selected column from the PDB dataframe
    df_concat_mdtraj_pdb = pd.concat([df_mdtraj_residues_string, extracted_columns], axis=1)
    
    
    
    #===== Create table of corresponding chainID --> chainName =====
    # create a dataframe with only chain ID and chain name
    df_chainID_chainName = df_concat_mdtraj_pdb[['chainID','chainName']]

    # remove duplicate rows
    df_chainID_chainName = df_chainID_chainName.drop_duplicates()
    
    
    
    #===== Generate dictionnary to convert chainID <--> chainName
    #----- chainID 2 chainName -----
    # generate the dictionnary from the dataframe 'df_chainID_chainName'
    dict_chainID_2_chainName = df_chainID_chainName.set_index('chainID')['chainName'].to_dict()
    
    # convert all chainID to int type
    dict_chainID_2_chainName = {int(key): value for key, value in dict_chainID_2_chainName.items()}
    
    
    #----- chainName 2 chainID -----
    # Initialize an empty reversed dictionary
    dict_chainName_2_chainID = {}

    # Iterate over the original dictionary
    for key, value in dict_chainID_2_chainName.items():
        # Append the value to the list associated with the key in the reversed dictionary
        dict_chainName_2_chainID.setdefault(value, []).append(int(key))
    
    
       
    #===== Create table of corresponding chainName --> chainID =====
    # Convert dictionary to DataFrame
    df_chainName_chainID = pd.DataFrame.from_dict(dict_chainName_2_chainID, orient='index')
    
    # rename all column to strat by 'chainID_'
    df_chainName_chainID.rename(columns={col: f'chainID_{col}' for col in df_chainName_chainID.columns}, inplace=True)
    
    # index (chainName) are set as column
    df_chainName_chainID.reset_index(inplace=True)
    df_chainName_chainID = df_chainName_chainID.rename(columns = {'index':'chainName'})
    

    
    #===== Save output files =====
    if write_outfile == True:
        # Get the output path with the original PDB name
        file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string

        # save the concatenate dataframe to a file
        df_concat_mdtraj_pdb.to_csv(f'{file_path}_mdtraj_pdb_correspondance.csv', index=False)
        
        # save the table of chainID corresponding to chainName
        df_chainID_chainName.to_csv(f'{file_path}_chainID_2_chainName.csv', index=False)
        
        # save the table of chainName corresponding to chainID 
        df_chainName_chainID.to_csv(f'{file_path}_chainName_2_chainID.csv', index=False)
        
        
    
    #===== Return convertion dictionnaries =====
    return dict_chainID_2_chainName, dict_chainName_2_chainID
    




#=====================================================
#===== Function to get sequence and secondary structure of each chain
#=====================================================
def get_sequence_secstruct(pdb_file_path, write_outfile=True):
    """
    DESCRIPTION
        Return the sequence and secondary structure of each chain in a PDB.

        The command returns two dictionaries:
        - chainID and sequence
        - chainID and secondary structure

        The command also writes 2 output files (as CSV):
        - correspondence between mdtraj toppology and PDB file
        - chainID and their corresponding sequence and secondary structure

        Please note that depending on the PDB file, MDTraj may use different
        IDs for the same chain name in the PDB file.
        
        
    USAGE
        chainID_sequence, chainID_secstructure = get_sequence_secstruct('5azz.pdb')
    
    
    ARGUMENT
        pdb_file_path    Path to the pdb file (string format)
        
        
    OPTIONAL ARGUMENT
        write_outfile    (True/False) Write the output files.   
    """
    #===== read the PDB file with MDTraj=====
    traj = md.load(pdb_file_path, top=pdb_file_path)


    #===== Get secondary structure =====
    # Compute DSSP vertion of MSTraj
    list_secondary_structure = md.compute_dssp(traj, simplified=True)[0]

    # replace all 'NA' to 'N'
    # it appear for << "residue" in the topology which isn't actually a protein residue. >>
    list_secondary_structure = list(map(lambda x: x.replace('NA', 'N'), list_secondary_structure))


    #===== Get residues parameters and COM, COM backbone, COM sidechain =====
    #----- initialise variables to store results -----
    list_name = []
    list_ID = []
    list_chainID = []
    list_COM = []
    list_COM_sidechain = []
    list_COM_backbone = []

    #----- loop over all residue in the topology -----
    for residue in traj.topology.residues:
        # append the list with the residue name
        list_name.append(residue.name)

        # append the list with the residue index
        list_ID.append(residue.index)

        # append the list with the residue index
        list_chainID.append(residue.chain.index)

        # calculate center of mass
        COM = md.compute_center_of_mass(traj, select=f"resid {residue.index}")[0]

        # calculate center of mass for sidechain only
        try:
            COM_sidechain = md.compute_center_of_mass(traj, select=f"resid {residue.index} and sidechain")[0]
        except:
            COM_sidechain = np.NaN

        # calculate center of mass for backbone only 
        try:
            COM_backbone = md.compute_center_of_mass(traj, select=f"resid {residue.index} and backbone")[0]
        except:
            COM_backbone = np.NaN

        # append the lists with the COMs
        list_COM.append(''.join(str(COM)))
        list_COM_sidechain.append(''.join(str(COM_sidechain)))
        list_COM_backbone.append(''.join(str(COM_backbone)))



    #===== Create a Pandas table =====
    df_topology = pd.DataFrame({'chainID':list_chainID,
                       'index':list_ID,
                       'name':list_name,
                       'secondary_structure':list_secondary_structure,
                       'COM':list_COM,
                       'COM_sidechain':list_COM_sidechain,
                       'COM_backbone':list_COM_backbone,     
                      })

    #===== 
    #----- Initialize variables -----
    list_sequence = []
    list_string_ss = []

    #
    for chainID in set(list_chainID):
        #..... select row containint the chainID in their column 'chainID' .....
        df_select_chainID = df_topology[(df_topology['chainID'].isin([chainID]))]

        #..... Get sequence of each chain .....
        # convert the column 'name' to list
        list_name_3_letter = df_select_chainID['name'].to_list()

        # convert the 3 letter code into 1 letter ode and stor it to another list
        list_name_1_letter = [seq1(name) for name in list_name_3_letter]

        # convert the 1 letter code list to a string
        sequence = ''.join(list_name_1_letter)

        #..... Get secondary structure sequence of the chain .....
        # convert the column 'secondary_structure' to list
        list_ss = df_select_chainID['secondary_structure'].to_list()

        # convert the secondary structure list to a string
        string_ss = ''.join(list_ss)

        #..... Append lists .....
        list_sequence.append(sequence)
        list_string_ss.append(string_ss)




    #===== Create pandas dataframe of chainID an their sequences =====
    df_chain_sequence_ss = pd.DataFrame({'chainID':list(set(list_chainID)),
                                         'sequence':list_sequence,
                                         'secondary_structure':list_string_ss,
                                        })




    #===== Save output files =====
    if write_outfile == True:
        # Get the output path with the original PDB name
        file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string

        # Save the corresponding 
        df_topology.to_csv(f'{file_path}_residues_secondaryStructure_and_COMs.csv', index=False)

        # Save the 
        df_chain_sequence_ss.to_csv(f'{file_path}_sequence_secondaryStructure.csv', index=False)





    #===== Generate output dictionnaries =====
    # generate the dictionnary chainID --> sequence
    dict_chainID_sequence = df_chain_sequence_ss.set_index('chainID')['sequence'].to_dict()

    # generate the dictionnary chainID --> secondary structure
    dict_chainID_ss = df_chain_sequence_ss.set_index('chainID')['secondary_structure'].to_dict()




    #===== Return result =====
    return dict_chainID_sequence, dict_chainID_ss





#=====================================================
#===== PDB Minimization
#=====================================================

#===== 1. Class to report minimisation steps =====
# It is not designed to be used by user
# Code modified from: https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/report_minimization.html

# The class can have any name but it must subclass MinimizationReporter.
class MyMinimizationReporter(MinimizationReporter):
    # within the class you can declare variables that persist throughout the minimization
    
    # you must override the report method and it must have this signature.
    def report(self, iteration, x, grad, args):
        '''
        the report method is called every iteration of the minimization.

        Args:
            iteration (int): The index of the current iteration. This refers
                             to the current call to the L-BFGS optimizer.
                             Each time the minimizer increases the restraint strength,
                             the iteration index is reset to 0.

            x (array-like): The current particle positions in flattened order:
                            the three coordinates of the first particle,
                            then the three coordinates of the second particle, etc.

            grad (array-like): The current gradient of the objective function
                               (potential energy plus restraint energy) with
                               respect to the particle coordinates, in flattened order.

            args (dict): Additional statistics described above about the current state of minimization.
                         In particular:
                         "system energy": the current potential energy of the system
                         "restraint energy": the energy of the harmonic restraints
                         "restraint strength": the force constant of the restraints (in kJ/mol/nm^2)
                         "max constraint error": the maximum relative error in the length of any constraint

        Returns:
            bool : Specify if minimization should be stopped.
        '''

        # Within the report method you write the code you want to be executed at
        # each iteration of the minimization.
        
        # The 'iteration', 'current energy', 'restraint energy', 'restraint strength' are wrote in an output file.            
        output_file = open('minimization_out.csv', "a")
        output_file.write(f"{iteration},{args['system energy']},{args['restraint energy']},{args['restraint strength']},{args['max constraint error']}\n")
        output_file.close()

        return False



#===== 2. Function to mnimize structure ===== 
def minimize_pdb(pdb_file_path, force_field='amber', max_iterations=100):
    """
    DESCRIPTION
        This command is a parser to minimize a structure using openMM.
        It generate a PDB file and an output file in CSV format containing
        the energy of the structure at each step of the minimization.

    ARGUMENTS
        pdb_file_path    Path of the pDB file.
        force_field      Force field 'amber' or 'charmm'.
                         Default value: 'amber'
        max_iterations   Maximum number of iteration. Too hight number can lead to structure deformation.
                         Default value: 100
    """
    #==== 
    #----- Clear existing output file and minimized pdb -----
    # remove the minimization output file if it exist
    if os.path.exists("minimization_out.csv"):
        os.remove("minimization_out.csv")
        
    # Get the output path with the original PDB name
    file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string
    # remove the minimized pdb if exist
    if os.path.exists(f'{file_path}_minimized.pdb'):
        os.remove(f'{file_path}_minimized.pdb')
    
    #----- create output file ----
    out = open('minimization_out.csv', "w")
    out.write("Iteration,System energy (kJ/mol),Harmonic restraints energy (kJ/mol),Restraint force constant (kJ/mol/nm^2),Constraint maximum relative error\n")
    out.close()
    
    
    
    #===== Initial papamerters =====
    # load the PDB file
    pdb = PDBFile(pdb_file_path)
    
    # get the forcefield
    if force_field == 'amber':
        forcefield = ForceField('amber14/protein.ff14SB.xml')
    elif force_field == 'charmm':
        forcefield = ForceField('charmm36.xml')
    else:
        raise ValueError('Forcefield value must be "amber" or "charmm"')
    
    
    #===== Set simulation parameter and minimize the structure =====
    # create the system
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)

    
    # specify integrator
    # paramters set are: temperature, friction coefficient, timestep
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    
    # generate a new simulation
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    
    # set minimization reporter, to export the log file
    reporter = MyMinimizationReporter()
    
    # perform the minimization
    simulation.minimizeEnergy(maxIterations=max_iterations, reporter=reporter)
    
    
    #===== Export minimized structure =====    
    # write the PDB file
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(f'{file_path}_minimized.pdb', 'w'))





#=====================================================
#===== Function to get SASA of resisues when bound or free in a complex
# Use in other functions. Not interesting to use as it by user.
#=====================================================
def get_SASA_residue_complex(pdb_file_path, chainID_receptor, chainID_ligand, ignore_hydrogen=False, dict_radii=None, probe_radius_angstrom=1.4, number_sphere_points=960):
    """
    DESCRIPTION
        Return the SASA of each residue in a protein complex when they are bonded or free.

    ARGUMENTS
        pdb_file_path       path of the PDB file
        chainID_receptor    list of MDTraj chainID of the receptor
        chainID_ligand      list of MDTraj chainID of the receptor

    OPTIONAL ARGUMENTS
        ignore_hydrogen     (True/False) Avoid hydrogen in SASA calculation.
                            Default vaule: False

        dict_radii          (False/dictionnary) Use customatom radii.
                            Default value: None.
                            More info: https://mdtraj.org/1.9.7/api/generated/mdtraj.shrake_rupley.html

        probe_radius_angstrom    Probe radius in angstrom.
                                 Default value: 1.4 

        number_sphere_points     number of points representing the surface of each atom
                                 Default value: 960                         

    """
    #===== read the PDB file with MDTraj=====
    traj = md.load(pdb_file_path, top=pdb_file_path)
    
    
    #===== Convert chainIDs to text =====
    # ensure chainID are string
    chainID_receptor = [str(i) for i in chainID_receptor]
    chainID_ligand = [str(i) for i in chainID_ligand]

    # join list of chainID as sigle string
    chainID_receptor = ' '.join(chainID_receptor)
    chainID_ligand = ' '.join(chainID_ligand)



    #===== Create traj of ligand, receptor =====
    # Remove or keep hydrogens atoms in the selections
    if ignore_hydrogen == False:
        if_hydrogens = ''
    elif ignore_hydrogen == True:
        if_hydrogens = 'and not element H'
    
    # Traj for the receptor
    select_receptor = traj.topology.select(f"chainid {chainID_receptor} and protein {if_hydrogens}") 
    traj_receptor = traj.atom_slice(select_receptor)

    # Traj for the ligand
    select_ligand = traj.topology.select(f"chainid {chainID_ligand} and protein {if_hydrogens}") 
    traj_ligand = traj.atom_slice(select_ligand)

    # Traj for the complex
    select_complex = traj.topology.select(f"chainid {chainID_ligand} {chainID_receptor} and protein {if_hydrogens}") 
    traj_complex = traj.atom_slice(select_complex)



    #===== dictionnary create dictionary to corellate CA position <--> residue index =====
    # dictionnary position --> original_resID
    dict_position_initial_indices = {" ".join(map(str,traj.xyz[0, atom_index])):traj.topology.atom(atom_index).residue.index for atom_index in traj.topology.select(f"chainid {chainID_ligand} {chainID_receptor} and name CA")}

    # dictionnary position --> ID in traj_complex --> position
    dict_complex_indices_positions = {" ".join(map(str,traj_complex.xyz[0, atom_index])):traj_complex.topology.atom(atom_index).residue.index for atom_index in traj_complex.topology.select(f"name CA")}

    # dictionnary position --> ID in traj_receptor --> position
    dict_receptor_indices_positions = {" ".join(map(str,traj_receptor.xyz[0, atom_index])):traj_receptor.topology.atom(atom_index).residue.index for atom_index in traj_receptor.topology.select(f"name CA")}

    # dictionnary position --> ID in traj_ligand --> position
    dict_ligand_indices_positions = {" ".join(map(str,traj_ligand.xyz[0, atom_index])):traj_ligand.topology.atom(atom_index).residue.index for atom_index in traj_ligand.topology.select(f"name CA")}



    #===== compute SASA =====
    # Convert probe_radius_angstrom in angstrom to nm
    probe_radius_angstrom = probe_radius_angstrom /10 

    # SASA of each residue in the traj_complex
    shrake_rupley_bound = md.shrake_rupley(traj_complex, mode='residue', get_mapping=True, change_radii=dict_radii, probe_radius=probe_radius_angstrom, n_sphere_points=number_sphere_points)
    SASA_bound = shrake_rupley_bound[0][0]
    SASA_bound_resID = list( set(shrake_rupley_bound[1]) )

    # SASA of each residue in the traj_receptor
    shrake_rupley_receptor = md.shrake_rupley(traj_receptor, mode='residue', get_mapping=True, change_radii=dict_radii, probe_radius=probe_radius_angstrom, n_sphere_points=number_sphere_points)
    SASA_receptor = shrake_rupley_receptor[0][0]
    SASA_receptor_resID = list( set(shrake_rupley_receptor[1]) )

    # SASA of each residue in the traj_ligand
    shrake_rupley_ligand = md.shrake_rupley(traj_ligand, mode='residue', get_mapping=True, change_radii=dict_radii, probe_radius=probe_radius_angstrom, n_sphere_points=number_sphere_points)
    SASA_ligand = shrake_rupley_ligand[0][0]
    SASA_ligand_resID = list( set(shrake_rupley_ligand[1]) )



    #===== Generate dictionaries to store resID --> SASA =====
    # for the complex
    dict_SASA_bound = {SASA_bound_resID[index]:SASA_bound[index] for index in range(len(SASA_bound_resID))}

    # for the receptor
    dict_SASA_receptor = {SASA_receptor_resID[index]:SASA_receptor[index] for index in range(len(SASA_receptor_resID))}

    # for the ligand
    dict_SASA_ligand = {SASA_ligand_resID[index]:SASA_ligand[index] for index in range(len(SASA_ligand_resID))}



    #===== Get SASA od  residue when in complex or free (not in complex) =====
    # create list to store resluts
    list_chainID = []
    list_index = []
    list_name = []
    list_SASA_bound = []
    list_SASA_free = []
    list_identity = []


    # loop over all CA position fin the complex
    for position in dict_complex_indices_positions:
        #----- Get original ID, name and chainID infos -----
        # use the CA position to get the original residue ID (from the traj)
        index_original_traj = dict_position_initial_indices[position]

        # get residue name from it's original index
        resname = traj.topology.residue(index_original_traj).name

        # get chainID from it's original index
        chainID = traj.topology.residue(index_original_traj).chain.index


        #----- Get SASA in complex -----
        # use the CA position to resID into the traj_complex
        index_complex_traj  = dict_complex_indices_positions[position]

        # Get the SASA of the residdue in complex
        SASA_bound = dict_SASA_bound[index_complex_traj]


        #----- Get SASA free -----  
        # use the CA position to resID into the traj_receptor (if in receptor) or traj_ligand (if in residue)
        # and get the SASA when not in complex (free)

        # try if in receptor
        try:
            index_receptor_traj = dict_receptor_indices_positions[position]
            SASA_free = dict_SASA_receptor[index_receptor_traj]
            identidy = "receptor"

        # try if in ligand
        except:    
            index_ligand_traj = dict_ligand_indices_positions[position]
            SASA_free = dict_SASA_ligand[index_ligand_traj]
            identidy = "ligand"



        #----- append lists with values -----
        list_index.append(index_original_traj)
        list_chainID.append(chainID)
        list_name.append(resname)
        list_SASA_bound.append(SASA_bound *100) # *100 to convert nm2 to A2
        list_SASA_free.append(SASA_free *100) # *100 to convert nm2 to A2
        list_identity.append(identidy)



    #===== Create Pandas dataframe to store info =====
    dataframe = pd.DataFrame({"identity": list_identity,
                              "chainID":list_chainID,
                              "resID":list_index,
                              "name":list_name,
                              "SASA_bound":list_SASA_bound,
                              "SASA_free":list_SASA_free,
                            })
    
    
    #===== Return results =====
    return dataframe






#=====================================================
#===== Function to get protein region of each amino acid in a protein complex
#=====================================================
def get_protein_region(pdb_file_path, chainID_receptor, chainID_ligand, write_outfile=True):
    """
    DESCRIPTION
        A command to return the protein region of all amino acids involved in a protein complex.
        It use diffrent maximum ASA values:
            - Thien et al. 2013 (https://doi.org/10.1371/journal.pone.0080635): Theoric
            - ibid. : Empiric
            - Miller et al. 1987 (https://doi.org/10.1016/0022-2836(87)90038-6)
            - Rose et al. 1985 (https://doi.org/10.1126/science.4023714)
            - Lins et al. 2003 (https://doi.org/10.1110/ps.0304803)
            - Samanta et al. 2002 (https://doi.org/10.1093/protein/15.8.659): Gly-X-Gly
            - ibid. : Ala-X-Ala
            - NACCESS software (http://www.ncbi.nlm.nih.gov/pubmed/994183, http://www.bioinf.manchester.ac.uk/naccess/)

    
    ARGUMENTS
        pdb_file_path       Path of the PDB file.
        chainID_receptor    list of MDTraj chain ID of the receptor.
        chainID_ligand      list of MDTraj chain ID of the ligand.
    
    OPTIONAL ARGUMENTS
        write_outfile       write the output CSV files.
    """
    
    #===== Define a dictionnary with all MaxASA value from bibliography =====
    dict_MaxASA = {
        # Thien et al. 2013 (https://doi.org/10.1371/journal.pone.0080635): Theoric
        'Thien_theoric': {'ALA': 129,
                          'ARG': 274,
                          'ASN': 195,
                          'ASP': 193,
                          'CYS': 167,
                          'GLU': 223,
                          'GLN': 225,
                          'GLY': 104,
                          'HIS': 224,
                          'ILE': 197,
                          'LEU': 201,
                          'LYS': 236,
                          'MET': 224,
                          'PHE': 240,
                          'PRO': 159,
                          'SER': 155,
                          'THR': 172,
                          'TRP': 285,
                          'TYR': 263,
                          'VAL': 174
                         },

        # Thien et al. 2013 (https://doi.org/10.1371/journal.pone.0080635): Empiric
        'Thien_empiric': {'ALA': 121,
                          'ARG': 265,
                          'ASN': 187,
                          'ASP': 187,
                          'CYS': 148,
                          'GLU': 214,
                          'GLN': 214,
                          'GLY': 97,
                          'HIS': 216,
                          'ILE': 195,
                          'LEU': 191,
                          'LYS': 230,
                          'MET': 203,
                          'PHE': 228,
                          'PRO': 154,
                          'SER': 143,
                          'THR': 163,
                          'TRP': 264,
                          'TYR': 255,
                          'VAL': 165
                         },

        # Miller et al. 1987 (https://doi.org/10.1016/0022-2836(87)90038-6): Empiric 
        'Miller': {'ALA': 113,
                   'ARG': 241,
                   'ASN': 158,
                   'ASP': 151,
                   'CYS': 140,
                   'GLU': 183,
                   'GLN': 189,
                   'GLY': 85,
                   'HIS': 194,
                   'ILE': 182,
                   'LEU': 180,
                   'LYS': 211,
                   'MET': 204,
                   'PHE': 218,
                   'PRO': 143,
                   'SER': 122,
                   'THR': 146,
                   'TRP': 259,
                   'TYR': 229,
                   'VAL': 160
                  },

        # Rose et al. 1985 (https://doi.org/10.1126/science.4023714): Empiric
        'Rose': {'ALA': 118.1,
                 'ARG': 256,
                 'ASN': 165.5,
                 'ASP': 158.7,
                 'CYS': 146.1,
                 'GLU': 186.2,
                 'GLN': 193.2,
                 'GLY': 88.1,
                 'HIS': 202.5,
                 'ILE': 181,
                 'LEU': 193.1,
                 'LYS': 225.8,
                 'MET': 203.4,
                 'PHE': 222.8,
                 'PRO': 146.8,
                 'SER': 129.8,
                 'THR': 152.5,
                 'TRP': 266.3,
                 'TYR': 236.8,
                 'VAL': 164.5
                },

        # Lins et al. 2003 (https://doi.org/10.1110/ps.0304803): Empiric
        'Lins': {'ALA': 111,
                 'ARG': 250,
                 'ASN': 166,
                 'ASP': 160,
                 'CYS': 157,
                 'GLN': 187,
                 'GLU': 194,
                 'GLY': 86,
                 'HIS': 191,
                 'ILE': 173,
                 'LEU': 179,
                 'LYS': 212,
                 'MET': 201,
                 'PHE': 208,
                 'PRO': 135,
                 'SER': 125,
                 'THR': 144,
                 'TRP': 249, 
                 'TYR': 227,
                 'VAL': 149
                },

        # Samanta et al. 2002 (https://doi.org/10.1093/protein/15.8.659): Empiric
        #     In the set 'Samanta_ala' the value for GLY come frome 'Samanta_gly' because their ise no data availble for Ala-Gly-Ala in the paper.
        'Samanta_gly': {'ALA': 116.4,
                        'ARG': 249.26,
                        'ASN': 168.87,
                        'ASP': 155.37,
                        'CYS': 141.48,
                        'GLN': 189.17,
                        'GLU': 187.16,
                        'GLY': 83.91,
                        'HIS': 198.51,
                        'ILE': 189.95,
                        'LEU': 197.99,
                        'LYS': 207.49,
                        'MET': 210.55,
                        'PHE': 223.29,
                        'PRO': 144.8,
                        'SER': 125.68,
                        'THR': 148.06,
                        'TRP': 265.42,
                        'TYR': 238.3,
                        'VAL': 162.24
                       },

        'Samanta_ala': {'ALA': 55.4,
                        'ARG': 190.24,
                        'ASN': 109.92,
                        'ASP': 97.8,
                        'CYS': 82.07,
                        'GLN': 129.68,
                        'GLU': 132.53,
                        'GLY': 83.91,
                        'HIS': 141.27,
                        'ILE': 130.71,
                        'LEU': 141.52,
                        'LYS': 147.99,
                        'MET': 150.39,
                        'PHE': 164.18,
                        'PRO': 106.44,
                        'SER': 69.08,
                        'THR': 88.62,
                        'TRP': 209.62,
                        'TYR': 180.03,
                        'VAL': 103.12
                       },

        # NACCESS software (http://www.ncbi.nlm.nih.gov/pubmed/994183, http://www.bioinf.manchester.ac.uk/naccess/)
        # Ala-X-Ala
        'NACCESS': {"ALA": 107.95,
                    "CYS": 134.28,
                    "ASP": 140.39,
                    "GLU": 172.25,
                    "PHE": 199.48,
                    "GLY": 80.10,
                    "HIS": 182.88,
                    "ILE": 175.12,
                    "LYS": 200.81,
                    "LEU": 178.63,
                    "MET": 194.15,
                    "ASN": 143.94,
                    "PRO": 136.13,
                    "GLN": 178.50,
                    "ARG": 238.76,
                    "SER": 116.50,
                    "THR": 139.27,
                    "VAL": 151.44,
                    "TRP": 249.36,
                    "TYR": 212.76
                    }    
    } # end of MaxASA dictionnary

    #===== =====
    dataframe = get_SASA_residue_complex(pdb_file_path, chainID_receptor, chainID_ligand)




    #===== Calculate the protein region using all method in dict_dict_MaxASA & expamd the dataframe =====
    for method in dict_MaxASA:
        #----- calculate rASA and drASA -----
        # calculate rSASA in complex or free using the MaxASA value for the method in dict_MaxASA
        dataframe[f'{method}_rASA_complex'] = dataframe[f'SASA_bound'] / dataframe['name'].map(dict_MaxASA[method])
        dataframe[f'{method}_rASA_free'] = dataframe[f'SASA_free'] / dataframe['name'].map(dict_MaxASA[method])

        # calculate drASA = rASA_free - rASA_complex
        dataframe[f'{method}_drASA'] = dataframe[f'{method}_rASA_free'] - dataframe[f'{method}_rASA_complex']

        # Initialize 'protein_region' column with default value
        dataframe[f'{method}_protein_region'] = np.nan

        #----- Identify protein region -----
        # surface & hydrated
        dataframe.loc[(dataframe[f'{method}_drASA'] == 0) & (40/100.0 < dataframe[f'{method}_rASA_complex']), f'{method}_protein_region'] = 'surface_hydrated'

        # surface
        dataframe.loc[(dataframe[f'{method}_drASA'] == 0) & (25/100.0 < dataframe[f'{method}_rASA_complex']) & (dataframe[f'{method}_rASA_complex'] <= 40/100.0), f'{method}_protein_region'] = 'surface'

        # interior
        dataframe.loc[(dataframe[f'{method}_drASA'] == 0) & (dataframe[f'{method}_rASA_complex'] < 25/100.0), f'{method}_protein_region'] = 'interior'

        # rim & NIS
        dataframe.loc[(dataframe[f'{method}_drASA'] != 0) & (dataframe[f'{method}_drASA'] <= 5/100.0) & (25/100.0 < dataframe[f'{method}_rASA_complex']), f'{method}_protein_region'] = 'rim_nis'    

        # rim & interaction
        dataframe.loc[(5/100.0 < dataframe[f'{method}_drASA']) & (25/100.0 < dataframe[f'{method}_rASA_complex']), f'{method}_protein_region'] = 'rim_interaction'

        # support
        dataframe.loc[(0 < dataframe[f'{method}_drASA']) & (dataframe[f'{method}_rASA_free'] < 25/100.0), f'{method}_protein_region'] = 'support'

        # core
        dataframe.loc[(0 < dataframe[f'{method}_drASA']) & (dataframe[f'{method}_rASA_complex'] < 25/100.0) & (25/100.0 < dataframe[f'{method}_rASA_free']), f'{method}_protein_region'] = 'core'



    #===== Create a dataframe with the total ASA of the complex, receptor, ligand and the area of the interface =====
    # calculate ASA of the complex
    ASA_total_complex = dataframe['SASA_bound'].sum()

    # calculate ASA of the receptor or ligand by filtering using the column 'identity'
    ASA_total_receptor = dataframe[dataframe['identity'] == 'receptor']['SASA_free'].sum()
    ASA_total_ligand = dataframe[dataframe['identity'] == 'ligand']['SASA_free'].sum()

    # calculate the interface
    interface = (ASA_total_receptor + ASA_total_ligand) - ASA_total_complex

    # create a dataframe with the results
    df_ASA_total_interface = pd.DataFrame({'identity':["complex", "receptor", "ligand", "interface"],
                                           'ASA':[ASA_total_complex, ASA_total_receptor, ASA_total_ligand, interface],
                                         })



    #===== Create dataframe to store protein region as sequence format =====
    # Get columns containing 'protein' in their title
    columns_protein_region = [col for col in dataframe.columns if 'protein_region' in col]

    # Create a new dataframe with selected column
    df_filtered_protein_region = dataframe[['identity', 'chainID'] + columns_protein_region]
    # Create a copy of the DataFrame
    df_filtered_protein_region = df_filtered_protein_region.copy()

    # Convert protein_region to their 1-letter code by replacing all occurrences
    df_filtered_protein_region.replace('surface_hydrated', 'H', inplace=True)
    df_filtered_protein_region.replace('surface', 'S', inplace=True)
    df_filtered_protein_region.replace('interior', 'I', inplace=True)
    df_filtered_protein_region.replace('rim_nis', 'N', inplace=True)
    df_filtered_protein_region.replace('rim_interaction', 'R', inplace=True)
    df_filtered_protein_region.replace('support', 'P', inplace=True)
    df_filtered_protein_region.replace('core', 'C', inplace=True)

    # create lists containing only chainID an their identityt (receptor or ligand)
    list_chainID = list(set(df_filtered_protein_region['chainID'].to_list()))
    list_chainID_identity = list(set(df_filtered_protein_region['identity'].to_list()))

    # initialise the new dataframe
    df_sequence_protein_region = pd.DataFrame({'chainID':list_chainID,
                                               'identity':list_chainID_identity,
                                              })

    # initialise a list to store all protein_region in sequence format
    list_sequence_protein_region = []

    # for each method (colum in the 'df_filtered_protein_region') convert get the
    for method in columns_protein_region:
        list_sequence = [] # (re)initialise the list for all method column

        # append the list_sequence for each chainID
        for ID in list_chainID:
            # selecti only row with the ID in their chainID column, then select only the column 'method', then convert it to list, finally joint the list
            sequece_protein_region = ''.join(df_filtered_protein_region[df_filtered_protein_region['chainID'] == ID][method].to_list())
            list_sequence.append(sequece_protein_region)

        # append the dataframe with a new colum corresponding to the method
        df_sequence_protein_region[method] = list_sequence

    
    
    #===== Save output files =====
    if write_outfile == True:
        # Get the output path with the original PDB name
        file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string
        
        # Save the dataframe with residue SASA, rASA, drASA, and protein region informations
        dataframe.to_csv(f'{file_path}_residues_protein_region.csv', index=False)
        
        # Save dataframe with area information of complex, free receptor, free ligand, and the interface receptor-ligand
        df_ASA_total_interface.to_csv(f'{file_path}_ASA_complex_receptor_ligand_and_interface.csv', index=False)
        
        # 
        df_sequence_protein_region.to_csv(f'{file_path}_chainID_protein_region_sequence_format.csv', index=False)
    
                            
        
    #===== Generate output dictionnaries =====
    # generate the dictionnary chainID --> sequence
    dict_ASA_total_interface = df_ASA_total_interface.set_index('identity')['ASA'].to_dict()

    # generate the dictionnary chainID --> secondary structure
    dict_sequence_protein_region = df_sequence_protein_region.set_index('chainID').to_dict()
                           
        
        
    #===== Return results =====
    return dict_ASA_total_interface, dict_sequence_protein_region    
    





#=====================================================
#===== Function to prepare the structure using PDB2PQR
#=====================================================
def pbd2pqr_parse(pdb_file_path, force_field='AMBER', ph=7.0, write_logfile=True):
    """
    DESCRITPTION

    ARGUMENTS

    OPTIONAL ARGUMENTS
    """
    # check force field
    if force_field not in ['AMBER','CHARMM']:
        raise ValueError('Force field must be AMBER or CHARMM')
        
    # Get the output path with the original PDB name
    file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string
    
    # Create output files as PQR and PDB format
    output_pqr_file = f"{file_path}.pqr"
    output_pdb_file = f"{file_path}_prq.pdb"
        
    # Run pdb2pqr as subprocess
    # Because << The [python] API is still changing and there is currently no guarantee that it will remain stable between minor releases. >>
    output = subprocess.run(['pdb2pqr', '--ff', force_field, '--ffout', force_field, '--titration-state-method', 'propka', '--with-ph', str(ph), '--drop-water', \
                    '--pdb-output', output_pdb_file, pdb_file_path, output_pqr_file], capture_output=True, text=True)
    
    
    # Save PDB2PQR output into log file
    if write_logfile == True:
        with open(f"{file_path}_pdb2pqr.log", "w") as file:
            file.write(output.stderr)
            file.write("="*100)
            file.write(output.stdout)
    





#=====================================================
#===== Import modules
#=====================================================