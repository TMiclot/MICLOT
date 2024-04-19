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



    #===== Generate output dictionnaries =====
    # generate the dictionnary chainID --> sequence
    dict_chainID_sequence = df_chain_sequence_ss.set_index('chainID')['sequence'].to_dict()

    # generate the dictionnary chainID --> secondary structure
    dict_chainID_ss = df_chain_sequence_ss.set_index('chainID')['secondary_structure'].to_dict()



    #===== Save output files =====
    if write_outfile == True:
        # Get the output path with the original PDB name
        file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string

        # Save the corresponding 
        df_topology.to_csv(f'{file_path}_residues_secondaryStructure_and_COMs.csv', index=False)

        # Save the 
        df_chain_sequence_ss.to_csv(f'{file_path}_sequence_secondaryStructure.csv', index=False)



    #===== Return result =====
    return dict_chainID_sequence, dict_chainID_ss





#=====================================================
#===== PDB Minimization
#=====================================================

#===== 1. Class to report minimisation steps =====
# It is not designed to be used by user
# Code come from: https://openmm.github.io/openmm-cookbook/dev/notebooks/cookbook/report_minimization.html

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
def minimize_pdb(pdb_file_path, ff='amber', max_iterations=100):
    """
    
    
    ! constraints method Hbonds lead to a multiple run of x iteration of minimisation in the output file. But conserve Hbond.
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
    if ff == 'amber':
        forcefield = ForceField('amber14/protein.ff14SB.xml')
    elif ff == 'charmm':
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
#===== Import modules
#=====================================================