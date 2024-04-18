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

import numpy as np
import mdtraj as md
import pandas as pd





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
        write_outfile    (True/False)
                         Write the correspondancs table between:
                         - correspondance between PDB file and MDTraj topology
                         - mdtraj chainID 2 PDB chainName
                         - PDB chainName 2 mdtraj chainID    
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
    topology = md.load(pdb_file, top=pdb_file).topology

    # Convert topology to dataframe
    df_mdtraj_residues, df_mdtraj_bonds = traj.topology.to_dataframe()

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
    df_chainName_chainID = pd.DataFrame.from_dict(reversed_dict, orient='index')
    
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
    