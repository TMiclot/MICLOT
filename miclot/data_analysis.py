#!/usr/bin/env python3

"""
This script is part of MICLOT ...
"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"


__all__ = ['plot', 'remove', 'clean_data']
# for readability classes call external functions, when it become too long to be read easly.
# all classes are at the end of this file.


# functions to generate clened data, final, tables, and analysis graphs, etc
# functions to delete generated csv (of cleaned data) and png files
# transformer le 'script_analysis.py' en une fonction de MICLOT ?


#=====================================================
#===== Import modules
#=====================================================
import os
import glob
from itertools import combinations
from collections import Counter
import logging
import ast
from tqdm import tqdm

from skspatial.objects import Vector
import seaborn as sns
import numpy as np

import pandas as pd
pd.options.mode.copy_on_write = True

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize, LogNorm







###################################################### General functions # dictionnary

#=====================================================
#===== Dictionnary to convert protein_region to their 1-letter code by replacing all occurrences
#=====================================================
dict_protein_region_1code = {
    'surface_hydrated': 'H', 
    'surface': 'S', 
    'interior': 'I', 
    'rim_nis': 'N', 
    'rim_interaction': 'R', 
    'support': 'P', 
    'core': 'C',    
    }




#=====================================================
#===== Function to convert string to numpy array
#=====================================================
def string_to_array(s):
    """
    Convert string like '[0.23 456.5 98.0]' to numpy array.
    """
    # Remove square brackets and split the string
    s = s.strip('[]')
    # Convert to numpy array
    array = np.array(s.split(), dtype=float)
    return array





#=====================================================
#===== Function to identify neighbor pairs based on neighbor residues
#=====================================================
def find_neighbor_couples(list_neighbor, list_couple):
    """
    list_neighbor = ["0_5", "3_4", "10_3", "1_9"]
    list_couple = ["0_1", "3_4", "10_12", "14_15", "3_5"]

    return dataframe like:
       pair_index  neighbor_pairs_index
    0  0_1         [3_5]
    1  3_4         [10_12, 3_5]
    2  10_12       [3_4, 3_5]
    3  14_15       NaN
    4  3_5         [0_1, 3_4, 10_12]
    """
    neighbor_map = {}
    
    # Parse the neighbors and create a dictionary to map each person to their neighbors
    for neighbors in list_neighbor:
        a, b = map(int, neighbors.split('_')) # convert "0_1" format to [0,1]
        if a not in neighbor_map:
            neighbor_map[a] = set()
        if b not in neighbor_map:
            neighbor_map[b] = set()
        neighbor_map[a].add(b)
        neighbor_map[b].add(a)
    
    data = []
    
    # Iterate through each couple and check for neighbors
    for couple_str in list_couple:
        a, b = map(int, couple_str.split('_')) # convert "0_1" format to [0,1]
        neighbor_couples = []
        for n_couple_str in list_couple:
            if n_couple_str == couple_str:
                continue
            x, y = map(int, n_couple_str.split('_'))
            if (x in neighbor_map.get(a, set()) or x in neighbor_map.get(b, set()) or
                y in neighbor_map.get(a, set()) or y in neighbor_map.get(b, set())):
                neighbor_couples.append(n_couple_str)
        if not neighbor_couples:
            neighbor_couples = np.nan
        data.append([couple_str, neighbor_couples])
    
    df = pd.DataFrame(data, columns=['pair_index', 'neighbor_pairs_index'])
    return df



#=====================================================
#===== Function to replace np.NaN by 'No' in a list (for pair neighbor)
#=====================================================
def map_list(lst, mapping_dict):
    """
    Return 'No' if the pair have no neighbors.
    Else return a list of item translated by the dicttoinnary mapping_dict
    """
    try:
        if pd.isna(lst):
            return "No"
    except:
        return [mapping_dict[item] for item in lst]




#=====================================================
#===== Function to to remove a file in a directory
#=====================================================
def rm(directory, name, file_format):
    """
    Fucntion to remove any files in a directory, base on a name and a file format.
    
    ARGUMENTS
        directory     path where to remole files
        name          complete name of motif contained in the name of files
        file_format   file format
    """
    # Find all 'format' files containing 'name' in their name 
    files_to_remove = glob.glob(os.path.join(directory, f'*{name}*.{file_format}'))
    
    # List to store name of removed files
    list_removed = []
    list_error = []
    # Remove each file
    for file_path in files_to_remove:
        try:
            os.remove(file_path)
            list_removed.append(file_path)
        except:
            list_error.append(file_path)
            continue
    return list_removed, list_error



#=====================================================
#===== Fuction to generate distance map using COMs
#=====================================================
def make_distance_map(directory, file_name_structure='residues_secondaryStructure_and_COMs', pdb_name=None, save=True):
    """
    DESCRIPTION
        Calculate all COMs distances and generate a CSV file with these values.

    ARGUMENTS
        directory     directory where CSV files are located
        
    OPTINAL ARGUMENTS
        file_name_structure    Name of file containing secondary structure and COMs information.
                               Default value: 'residues_secondaryStructure_and_COMs'
                                    
        pdb_name    structure name to use in the clened file: {pdb_name}_distance_map.csv
                    Default value: directory name
                    
        save        save cleaned file as CSV in the directory
                    Default value: True
    """
    #===== Read structure file as dataframe =====
    file_structure = glob.glob(os.path.join(directory, f'*{file_name_structure}*.csv'))[0]
    df_structure = pd.read_csv(file_structure)
    
    
    #===== Create a list to store the distance data ====
    list_distance = []

    
    #===== Calculates COMs distances and append the list_distance =====
    for (index_residue_1, com_residue_1, com_sidechain_1, com_backbone_1), (index_residue_2, com_residue_2, com_sidechain_2, com_backbone_2) in combinations(zip(df_structure['index'], df_structure['COM'], df_structure["COM_sidechain"], df_structure["COM_backbone"]), 2):
        
        #----- Convert COMS position from string to numpy array -----
        com_residue_1   = string_to_array(com_residue_1)
        com_residue_2   = string_to_array(com_residue_2)
        com_sidechain_1 = string_to_array(com_sidechain_1)
        com_sidechain_2 = string_to_array(com_sidechain_2)
        com_backbone_1  = string_to_array(com_backbone_1)
        com_backbone_2  = string_to_array(com_backbone_2)
        
        #----- Calulate distances: COM-COM, COMbackbone-COMbackbone, COMsidechain-COMsidechain -----
        distance_com = Vector.from_points(com_residue_1, com_residue_2).norm() *10 #*10 to convert nm to angstrom
        distance_com_backbone  = Vector.from_points(com_backbone_1, com_backbone_2).norm() *10 #*10 to convert nm to angstrom
        distance_com_sidechain = Vector.from_points(com_sidechain_1, com_sidechain_2).norm() *10 #*10 to convert nm to angstrom
        
        #----- Calulate distances: COMbackbone-COMsidechain -----
        distance_com_backbone_residue_1_sidechain_residue_2 = Vector.from_points(com_backbone_1, com_sidechain_2).norm() *10 #*10 to convert nm to angstrom
        distance_com_backbone_residue_2_sidechain_residue_1 = Vector.from_points(com_backbone_2, com_sidechain_1).norm() *10 #*10 to convert nm to angstrom
        
        #----- Calulate distances: COM-COMbackbone, COM-COMsidechain -----
        distance_com_residue_1_backbone_residue_2  = Vector.from_points(com_residue_1, com_backbone_2).norm()  *10 #*10 to convert nm to angstrom
        distance_com_residue_1_sidechain_residue_2 = Vector.from_points(com_residue_1, com_sidechain_2).norm() *10 #*10 to convert nm to angstrom
        distance_com_residue_2_backbone_residue_1  = Vector.from_points(com_residue_2, com_backbone_1).norm()  *10 #*10 to convert nm to angstrom
        distance_com_residue_2_sidechain_residue_1 = Vector.from_points(com_residue_2, com_sidechain_1).norm() *10 #*10 to convert nm to angstrom
        
        #----- Append list of distance with new informations -----
        list_distance.append({'residue_1_index': index_residue_1,
                              'residue_2_index': index_residue_2,
                              'distance_COM_1_COM_2': distance_com,
                              'distance_COM_backbone_1_COM_backbone_2':   distance_com_backbone,
                              'distance_COM_sidechain_1_COM_sidechain_2': distance_com_sidechain,
                              'distance_COM_backbone_1_COM_scidechain_2': distance_com_backbone_residue_1_sidechain_residue_2,
                              'distance_COM_backbone_2_COM_scidechain_1': distance_com_backbone_residue_2_sidechain_residue_1,
                              'distance_COM_1_COM_backbone_2':   distance_com_residue_1_backbone_residue_2,
                              'distance_COM_1_COM_scidechain_2': distance_com_residue_1_sidechain_residue_2,
                              'distance_COM_2_COM_backbone_1':   distance_com_residue_2_backbone_residue_1,
                              'distance_COM_2_COM_scidechain_1': distance_com_residue_2_sidechain_residue_1,
                             })

        
    #===== Create a new DataFrame from the distance data =====
    df_distances = pd.DataFrame(list_distance)

    
    #===== Create pair index code (ex: 5_16) =====
    df_distances["pair_index"] = df_distances.apply(
            lambda row: '_'.join(map(str, sorted([row['residue_1_index'].astype(int), row['residue_2_index'].astype(int)]))),
            axis=1
        )
    
    
    #===== Check if residues are consecutive =====
    # consecutive residue are: i, i+1
    df_distances['consecutive'] = (abs(df_distances['residue_1_index'] - df_distances['residue_2_index']) == 1).astype(int)

    
    #===== save to CSV =====
    if save == True:

        if pdb_name == None:
            pdb_name = directory.split('/')[-1]
            
        path_clean_distance_map = f'{directory}/{pdb_name}_distance_map.csv'
        df_distances.to_csv(path_clean_distance_map, index=False)
    
    #===== return dataframe =====
    return df_distances





###################################################### Cleaning & Concatenate
#=====================================================
#===== Function to generate neigbor residues
#=====================================================
def clean_neighbor_residues(directory, file_name_clean_structure='clean_structure', \
                    file_name_distance_map='distance_map', neighbor_cutoff_distance=8, distance_columns=None, pdb_name=None, save=True):
    """
    DESCRIPITON
        Use the distance map to identify neighbor residues based on their COMs distances.
        
    ARGUMENTS
        directory     directory where CSV files are located
        
    OPTIONAL ARGUMENTS
        file_name_structure    Name of file containing secondary structure and COMs information.
                               Default value: 'residues_secondaryStructure_and_COMs'
                               
        neighbor_cutoff_distance    Cutoff distance used to identify neighbor residues.
                                    If any value in the list of columns distance_columns' is below or equal
                                    to this distance, the residues are neighbors.
                                    Default value: 8 angstrom
        
        distance_columns    list of columns in which to search for the cutoff distances
                            Default value: None (search in all distances columns)
        
        pdb_name    structure name to use in the clened file: {pdb_name}_distance_map.csv
                    Default value: directory name
                    
        save        save cleaned file as CSV in the directory
                    Default value: True 
    """
    #===== Get structure information =====
    #----- Read CSV files as pandas dataframe -----
    file_structure = glob.glob(os.path.join(directory, f'*{file_name_clean_structure}*.csv'))[0]
    df_structure = pd.read_csv(file_structure)
    
    #----- Conversion dictionnaries -----
    dict_residue_index2codeComplete   = df_structure.set_index('index')['code_complete'].to_dict()
    dict_residue_index2codeSimplified = df_structure.set_index('index')['code_name_secondary_structure'].to_dict()
    dict_residue_index2name           = df_structure.set_index('index')['name'].to_dict()
    dict_residue_index2ss             = df_structure.set_index('index')['secondary_structure'].to_dict()
    dict_residue_index2pr             = df_structure.set_index('index')['protein_region'].to_dict()
    dict_residue_index2identity       = df_structure.set_index('index')['identity'].to_dict()

    
    #===== Read distance map file, or create it =====
    try:
        file_distance_map = glob.glob(os.path.join(directory, f'*{file_name_distance_map}*.csv'))[0]
        df_distances = pd.read_csv(file_distance_map)
    except:
        df_distances = make_distance_map()

    
    #===== Create neighbor df based on the neighbor_cutoff_distance value =====
    if distance_columns == None:
        distance_columns = df_distances.filter(like='distance').columns # use all column containing 'distance' in their name

    # select row containing <= neighbor_cutoff_distance in any of their 'distance_columns'
    condition = df_distances[distance_columns].le(neighbor_cutoff_distance).any(axis=1) 

    # Create the neighbor residue dataframe
    df_neighbor_residues = df_distances[condition].copy()

    
    #===== Residues information: codes, names, etc =====
    # Codes
    df_neighbor_residues["residue_1_codeComplete"]   = df_neighbor_residues['residue_1_index'].map(dict_residue_index2codeComplete)
    df_neighbor_residues["residue_2_codeComplete"]   = df_neighbor_residues['residue_2_index'].map(dict_residue_index2codeComplete)
    df_neighbor_residues["residue_1_codeSimplified"] = df_neighbor_residues['residue_1_index'].map(dict_residue_index2codeSimplified)
    df_neighbor_residues["residue_2_codeSimplified"] = df_neighbor_residues['residue_2_index'].map(dict_residue_index2codeSimplified)

    # names
    df_neighbor_residues["residue_1_name"] = df_neighbor_residues['residue_1_index'].map(dict_residue_index2name)
    df_neighbor_residues["residue_2_name"] = df_neighbor_residues['residue_2_index'].map(dict_residue_index2name)

    # secondary structures
    df_neighbor_residues["residue_1_secondary_structure"] = df_neighbor_residues['residue_1_index'].map(dict_residue_index2ss)
    df_neighbor_residues["residue_2_secondary_structure"] = df_neighbor_residues['residue_2_index'].map(dict_residue_index2ss)

    # protein region
    df_neighbor_residues["residue_1_protein_region"] = df_neighbor_residues['residue_1_index'].map(dict_residue_index2pr)
    df_neighbor_residues["residue_2_protein_region"] = df_neighbor_residues['residue_2_index'].map(dict_residue_index2pr)

    # residues identity (receptor/ligand)
    df_neighbor_residues["residue_1_identity"] = df_neighbor_residues['residue_1_index'].map(dict_residue_index2identity)
    df_neighbor_residues["residue_2_identity"] = df_neighbor_residues['residue_2_index'].map(dict_residue_index2identity)

    
    #===== Replace all np.NaN values by string 'NaN' =====
    df_neighbor_residues = df_neighbor_residues.fillna("NaN")

    
    #===== Generate pair codes =====
    # Create the "neighbor_code_full" column with sorted combined information
    df_neighbor_residues["neighbor_code_full"] = df_neighbor_residues.apply(
            lambda row: '-'.join(map(str, sorted([row['residue_1_codeComplete'], row['residue_2_codeComplete']]))),
            axis=1
        )

    # Create the "neighbor_code_full" column with sorted combined information
    df_neighbor_residues["neighbor_code_name_secondary_structure"] = df_neighbor_residues.apply(
            lambda row: '-'.join(map(str, sorted([row['residue_1_codeSimplified'], row['residue_2_codeSimplified']]))),
            axis=1
        )

    # Create code for NAME column with sorted combined information
    df_neighbor_residues["neighbor_code_name"] = df_neighbor_residues.apply(
            lambda row: '-'.join(map(str, sorted([row['residue_1_name'], row['residue_2_name']]))),
            axis=1
        )

    # Create code for SS column with sorted combined information
    df_neighbor_residues["neighbor_code_secondary_structure"] = df_neighbor_residues.apply(
            lambda row: '-'.join(map(str, sorted([row['residue_1_secondary_structure'], row['residue_2_secondary_structure']]))),
            axis=1
        )

    # Create code for PR column with sorted combined information
    df_neighbor_residues["neighbor_code_protein_region"] = df_neighbor_residues.apply(
            lambda row: '-'.join(map(str, sorted([row['residue_1_protein_region'], row['residue_2_protein_region']]))),
            axis=1
        )


    #===== Create the identity of the pair =====
    df_neighbor_residues["neighbor_identity"] = df_neighbor_residues.apply(
            lambda row: '_'.join(map(str, sorted([row['residue_1_identity'], row['residue_2_identity']]))),
            axis=1
        )


    #===== Remove unwanted columns =====
    columns_to_remove = ["residue_1_index", "residue_2_index",
             "residue_1_name", "residue_2_name",
             "residue_1_codeComplete", "residue_2_codeComplete",
             "residue_1_codeSimplified", "residue_2_codeSimplified",
             "residue_1_secondary_structure", "residue_2_secondary_structure",
             "residue_1_protein_region", "residue_2_protein_region",
             "residue_1_identity", "residue_2_identity",
            ]

    df_neighbor_residues = df_neighbor_residues.drop(columns_to_remove, axis=1)
    df_neighbor_residues.reset_index(inplace=True, drop=True)

    
    #===== save to CSV =====
    if save == True:

        if pdb_name == None:
            pdb_name = directory.split('/')[-1]
            
        path_file = f'{directory}/{pdb_name}_clean_neighbor_residues.csv'
        df_neighbor_residues.to_csv(path_file, index=False)

    
    #===== return dataframe =====
    return df_neighbor_residues





#=====================================================
#===== Function to clean structure information
#=====================================================
def clean_structure(directory, file_name_structure='residues_secondaryStructure_and_COMs', \
                    file_name_protein_region='residues_protein_region', pdb_name=None, save=True):
    """
    DESCRIPTION
        Use CSV files: *residues_secondaryStructure_and_COMs.csv
                       *residues_protein_region.csv    (if it don't exist, protein_region and identity will be replaced by 'No')
        
        Clean structure information and write a table containing:
            chainID, index, name, secondary structure,
            protein region (for complex), identity (for complex),
            code complete, code simplified,
            distances: COM_COMsidechain, COM_COMbackbone, COMsidechain_COMbackbone
    
    ARGUMENTS
        directory     directory where CSV files are located
        
    OPTINAL ARGUMENTS
        file_name_structure         Name of file containing secondary structure and COMs information.
                                    Default value: 'residues_secondaryStructure_and_COMs'
                               
        file_name_protein_region    Name of file containing protein region and identity information.
                                    Default value: 'residues_protein_region'
                                    
        pdb_name    structure name to use in the clened file: {pdb_name}_clean_structure.csv
                    Default value: directory name
                    
        save        save cleaned file as CSV in the directory
                    Default value: True
    """    
    
    #===== Read structure file as dataframe =====
    file_structure = glob.glob(os.path.join(directory, f'*{file_name_structure}*.csv'))[0]
    df_structure = pd.read_csv(file_structure) # will be use as final dataframe
    
    
    #===== Read protein region file as dataframe, if any =====
    try:
        file_protein_region = glob.glob(os.path.join(directory, f'*{file_name_protein_region}*.csv'))[0]
        df_protein_region = pd.read_csv(file_protein_region)
    
    except:
        file_protein_region = None
        
        
    #===== Add protein region and identity information, if any =====
    if file_protein_region != None:
        # Select all columns containing "_protein_region"
        columns_protein_region = df_protein_region.filter(like='_protein_region')
        
        # Find the most frequent text in each row across the "_protein_region" columns
        df_protein_region['protein_region'] = columns_protein_region.mode(axis=1)[0]  
        
        # Get 1 letter code of protein region and identity
        df_structure['protein_region'] = df_protein_region['protein_region'].map(dict_protein_region_1code)
        df_structure['identity'] = df_protein_region['identity']
    
    else:
        df_structure['protein_region'] = 'No'
        df_structure['identity'] = 'No'
        
        
    #===== Generate residue codes =====
    df_structure["code_complete"] = df_structure["name"] + "_" + df_structure["secondary_structure"] + "_" + df_structure['protein_region']
    df_structure["code_name_secondary_structure"] = df_structure["name"] + "_" + df_structure["secondary_structure"]
    
    
    #===== Compute distances =====
    # Convert COM column string values as array values
    df_structure["COM"] = df_structure["COM"].map(string_to_array)
    df_structure["COM_sidechain"] = df_structure["COM_sidechain"].map(string_to_array)
    df_structure["COM_backbone"] = df_structure["COM_backbone"].map(string_to_array)
    
    # Compute distances: lenght of the vector between the COM points
    df_structure["distance_COM_backbone_COM_sidechain"] = df_structure.apply(
        lambda row: Vector.from_points(row["COM_backbone"], row["COM_sidechain"]).norm() *10, #*10 to convert nm to angstrom,
        axis=1
    )

    df_structure["distance_COM_COM_sidechain"] = df_structure.apply(
            lambda row: Vector.from_points(row["COM"], row["COM_sidechain"]).norm() *10, #*10 to convert nm to angstrom,
            axis=1
        )

    df_structure["distance_COM_COM_backbone"] = df_structure.apply(
            lambda row: Vector.from_points(row["COM"], row["COM_backbone"]).norm() *10, #*10 to convert nm to angstrom,
            axis=1
        )
    
    
    #===== Remove unwanted columns =====
    columns_to_remove = ["COM",
                         "COM_backbone",
                         "COM_sidechain"
                        ]

    df_structure = df_structure.drop(columns_to_remove, axis=1)
   
    
    #===== Save as CSV =====
    if save == True:

        if pdb_name == None:
            pdb_name = directory.split('/')[-1]

        path_clean_df = f'{directory}/{pdb_name}_clean_structure.csv'
        df_structure.to_csv(path_clean_df, index=False)
    
    #===== Return final dataframe =====
    return df_structure
    






#=====================================================
#===== Function to concatenate CSV files
#=====================================================
def concatenate_csv(directory, file_name, final_file_name=None, logfile_name=None, use_tqdm=False):
    """
    DESCRIPTION
        Concatenade CSV files.
    
    ARGUMENTS
        directory    directory containing all subdirectory where to find CSV files
        
        file_name    complete file name or motif in the name of files to concatenate
        
    OPTIONAL ARGUMENTS
        final_file_name    name of CSV file corresponding to concatenation
        
        logfile_name       name of the log file
        
        use_tqdm           Use or not TQDM. (True/False)
                           Default value: 
    
    """
    #===== Configure log file =====
    if logfile_name is None:
        logfile_name = "concatenate"
        
    # Set log file path
    logfile_path = os.path.join(directory, f"{logfile_name}.log")

    # Create a unique logger name based on the logfile name
    logger_name = f"logger_{logfile_name}"
    logger = logging.getLogger(logger_name)
    
    # If the logger already has handlers, remove them to avoid duplicate logs
    if logger.hasHandlers():
        logger.handlers.clear()
    
    # Create file handler
    file_handler = logging.FileHandler(logfile_path, encoding='utf-8')
    
    # Create formatter and set it for the handler
    formatter = logging.Formatter('%(asctime)s - %(levelname)s :: %(message)s')
    file_handler.setFormatter(formatter)
    
    # Add handler to the logger
    logger.addHandler(file_handler)
    
    # Set the logging level
    logger.setLevel(logging.INFO)
    
    # Append log file
    logger.info(">>> START the CSV files concatenation <<<")
    
    
    #===== Remove '.csv' in 'file_name' variable if any =====
    if file_name.endswith('.csv'):
        file_name = file_name.split('.')[0]
    
    
    #===== Set final_file_name =====
    if final_file_name == None:
        final_file_name = f'FINAL_{file_name}.csv'
    
    final_file_name = os.path.join(directory, final_file_name)
    
    
    #===== Get list of all sub-directories =====
    dict_subdir = {root:subdir for root, subdir, _ in os.walk(directory)}
    list_subdir = [i for i in dict_subdir[directory] if not i.startswith('.')] # remove all hiden subdirectory
    
    #----- Setup progress for log file -----
    total_dirs = len(list_subdir)
    processed_dirs = 0
    
    
    #===== Concatenate =====
    if use_tqdm == True:
        list_subdir = tqdm(list_subdir)
    
    for subdir in list_subdir:
        #----- Update progress in log file -----
        processed_dirs += 1
        progress_percent = float((processed_dirs/total_dirs) *100.0)
        
        # Write progress in the log file
        logger.info(f"> Working on: {subdir} | Progress: {processed_dirs}/{total_dirs} ({progress_percent} %)")
        
        #----- generate the complete path of the subdir -----
        subdir_path = os.path.join(directory, subdir)

        #----- set pdb_name -----
        pdb_name = subdir

        #----- Get CSV file containing file_name -----
        # Get CSV file path
        csv_file = glob.glob(os.path.join(subdir_path, f'*{file_name}*.csv'))[0]
        
        # read file as dataframe
        df_table = pd.read_csv(csv_file)
        
        # transpose table for ASA
        if "ASA_complex_receptor_ligand_and_interface" in csv_file:
            df_table = df_table.set_index('identity').T.reset_index(drop=True)
        
        #----- add PDB name -----
        df_table["PDB"] = pdb_name
        
        #----- read final file or create it -----
        try:
            df_final_table = pd.read_csv(final_file_name)
            
            # concatenate
            df_final_table = pd.concat([df_table, df_final_table], axis=0)
            logger.info(f"Concatenate: {csv_file} with {final_file_name}")
            
        except FileNotFoundError:
            # if final file d'ont exist, use df_table
            df_final_table = df_table
            
        except EmptyDataError:
            logger.warning(f"Unable to read CSV files: {final_file_name}")
            continue
            
        #----- save CSV -----
        df_final_table.to_csv(final_file_name, index=None)
        logger.info(f"Save CSV file: {final_file_name}") 
       
    #===== End of normal execution =====
    logger.info("> Stop normaly")
    logger.info(">>> END of CSV file concatenation <<<")
    logging.shutdown()

    return concatenate_csv




###################################################### Plots

#=====================================================
#===== Function to plot minimization energies
#=====================================================
def plot_minimization(path, name=None, fig_size=[15, 4], linewidth=2, save=False, dpi=300):
    """
    DESCRIPTION
        Function used to plot minimisation energies from the CSV output file.

        
    ARGUMENTS
        path      directory where to find the '*_minimization_log.csv' file.
                  Alternatively it can be also the path of the file.

                  
    OPTIONAL ARGUMENTS
        name    By defaut the directory name where the CSV file is located.
                    Example 'my super PDB'
                    Default value: None
                  
        fig_size    list contiaing figure size [x_size, y_sizer]
                    Default value: [15,4]
        
        linewidth    line width of the graph.
                     Default value: 2
                   
        save    save the graph in the same directory as the CSV file. (True/False)
                      Default value: False
                      
        dpi    Dots per Inch. When save the graph.
               Default value: 300 
    """
    #===== Dictionnanry to store title names =====
    dict_title = {'system_energy':'System energy',
                  'harmonic_restraints_energy':'Harmonic restraints energy',
                  'restraint_force_constant':'Restraint force constant',
                  'constraint_maximum_relative_error':'Constraint maximum relative error',
                 }
    
    
    #===== Get the outputed CSV file of the minimisation =====
    if path.endswith(".csv"):
        minimization_log = path
        
    else:
        #Automaticaly search for the energy log file in path
        file_pattern = os.path.join(path, f'*_minimization_log.csv')
        minimization_log = glob.glob(file_pattern)[0]
    
    # Sample DataFrame creation for illustration purposes
    df = pd.read_csv(f'{minimization_log}')
    
    
    #==== Get PDB name ====
    if name == None:
        if path.endswith(".csv"): # if file path is provided with file name
            name = path.split('/')[-2]
        else: # if directory path is provided (no file name)
            name = path.split('/')[-1]
        
    
    #===== Create the plot =====    
    # Create 4 subgraph for energies
    fig, ax = plt.subplot_mosaic([list(dict_title)], layout='constrained', figsize=(fig_size[0],fig_size[1]))
    
    for column in dict_title:
        ax[f'{column}'].plot(df[f'{column}'], linewidth=linewidth)

        # Set title and labels
        ax[f'{column}'].set_title(dict_title[column])
        ax[f'{column}'].set_xlabel('Iteration')
        ax[f'{column}'].set_ylabel('Energy (kJ/mol)')

        
    #===== Set graph title =====
    fig.suptitle(f'Energy variation during minimization - {name}', fontsize=20)

    #===== Save the plot as a high-quality PNG file =====
    if save == True:
        plot_name = minimization_log.split('.')[0]
        plt.savefig(f'{plot_name}_plot_minimization_energies.png', dpi=dpi)
    
    #===== return the plot =====
    return plt



###################################################### Classes

#=====================================================
#===== Class to remove files
#=====================================================
class remove():
    def __init__(self):
        """
        Contain functions to remove generated cleaned data and plots.
        """
        pass

    #===== Remove any file generated by cleaning data =====
    def files(*arg):
        """
        Fucntion to remove any files in a directory, base on a name and a file format.
        
        ARGUMENTS
            directory     path where to remole files
            name          complete name of motif contained in the name of files
            file_format   file format
        """
        list_removed, list_error = rm(*arg)
        return list_removed, list_error

    #===== Remove CSV file generated by cleaning data =====
    def clean_data(*arg):
        """
        In a directory, remove all CSV files containing 'clean' in their name.
        """
        list_removed, list_error = rm(*arg, name='clean', file_format='csv')
        return list_removed, list_error


    #===== Remove PNG file generated by cleaning data =====
    def plots(*arg):
        """
        In a directory, remove all PNG files containing 'plot' in their name.
        """
        list_removed, list_error = rm(*arg, name='plot', file_format='png')
        return list_removed, list_error
    

    #===== Remove PNG file generated by cleaning data =====
    def pickles(*arg):
        """
        In a directory, remove all compressed pickle files containing 'class' in their name.
        """
        list_removed, list_error = rm(*arg, name='class', file_format='pkl.gz')
        return list_removed, list_error





#=====================================================
#===== Class to plot data
#=====================================================
class plot():
    def __init__(self):
        """
        Contain functions to easely plot data.
        """
        pass

    #===== Plot energy from minimization output ===== 
    def minimization(*arg,**kwargs):
        """
        DESCRIPTION
            Function used to plot minimisation energies from the CSV output file.
            
        ARGUMENTS
            path    directory where to find the '*_minimization_log.csv' file.
                    Alternatively it can be also the path of the file.
          
        OPTIONAL ARGUMENTS
            name    By defaut the directory name where the CSV file is located.
                        Example 'my super PDB'
                        Default value: None
                    
            fig_size    list contiaing figure size [x_size, y_sizer]
                        Default value: [15,4]
            
            linewidth    line width of the graph.
                        Default value: 2
                    
            save    save the graph in the same directory as the CSV file. (True/False)
                        Default value: False
                        
            dpi    Dots per Inch. When save the graph.
                Default value: 300 
        """
        return plot_minimization(*arg,**kwargs)





#=====================================================
#===== Class to concatenate (cleaned) CSV files
#=====================================================
class concatenate:
    def __init__():
        """
        Functions used to concatenate CSV files
        """
        pass
    
    #===== Concatenate CSV files over subdirectory ===== 
    def csv_file(*arg):
        """
        DESCRIPTION
            Concatenade CSV files.
        
        ARGUMENTS
            directory    directory containing all subdirectory where to find CSV files
            
            file_name    complete file name or motif in the name of files to concatenate
            
        OPTIONAL ARGUMENTS
            final_file_name    name of CSV file corresponding to concatenation
            
            logfile_name       name of the log file
            
            use_tqdm           Use or not TQDM. (True/False)
                            Default value: 
        
        """
        concatenate_csv(*arg)
        

    #===== Concatenate files with ASA information =====
    def ASA(*arg, **kwargs):
        """
        Concatenate '*ASA_complex_receptor_ligand_and_interface.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'ASA_complex_receptor_ligand_and_interface.csv', 'final_ASA_complexes.csv', logfile_name='concatenate_ASA', **kwargs)
        

    #===== Concatenate files with interaction_table information =====
    def interaction_table(*arg, **kwargs):
        """
        Concatenate '*clean_interactions_table.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'clean_interactions_table.csv', 'final_interactions_table.csv', logfile_name='concatenate_interactions_table', **kwargs)
        

    #===== Concatenate files with neighbor_residues information =====
    def neighbor_residues(*arg, **kwargs):
        """
        Concatenate '*clean_neighbor_residues.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'clean_neighbor_residues.csv', 'final_neighbor_residues.csv', logfile_name='concatenate_neighbor_residues', **kwargs)
           

    #===== Concatenate files with neighbor_pairs information =====
    def neighbor_pairs(*arg, **kwargs):
        """
        Concatenate '*clean_neighbor_pairs.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'clean_neighbor_pairs.csv', 'final_neighbor_pairs.csv', logfile_name='concatenate_neighbor_pairs', **kwargs)
        

    #===== Concatenate files with structure information =====
    def structure(*arg, **kwargs):
        """
        Concatenate '*clean_structure.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'clean_structure.csv', 'final_structure.csv', logfile_name='concatenate_', **kwargs)





#=====================================================
#===== Class to clean data
#=====================================================
class cleaning:
    def __init__():
        """
        Functions used to clean data
        """
        pass
    
    def all(*arg):
        """
        Perform all cleaning procedure in a directory. Save everything as CSV and don't return dataframes.
        """
        clean_structure(*arg)

    #===== Clean structure info ===== 
    def structure(*arg,**kwargs):
        """
        DESCRIPTION
            Use CSV files: *residues_secondaryStructure_and_COMs.csv
                        *residues_protein_region.csv    (if it don't exist, protein_region and identity will be replaced by 'No')
            
            Clean structure information and write a table containing:
                chainID, index, name, secondary structure,
                protein region (for complex), identity (for complex),
                code complete, code simplified,
                distances: COM_COMsidechain, COM_COMbackbone, COMsidechain_COMbackbone
        
        ARGUMENTS
            directory     directory where CSV files are located
            
        OPTINAL ARGUMENTS
            file_name_structure         Name of file containing secondary structure and COMs information.
                                        Default value: 'residues_secondaryStructure_and_COMs'
                                
            file_name_protein_region    Name of file containing protein region and identity information.
                                        Default value: 'residues_protein_region'
                                        
            pdb_name    structure name to use in the clened file: {pdb_name}_clean_structure.csv
                        Default value: directory name
                        
            save        save cleaned file as CSV in the directory
                        Default value: True
        """
        return clean_structure(*arg,**kwargs)


    #===== Generate clean neighbor residue using distance map ===== 
    def neighbor_residues(*arg,**kwargs):
        """
        DESCRIPITON
            Use the distance map to identify neighbor residues based on their COMs distances.
            
        ARGUMENTS
            directory     directory where CSV files are located
            
        OPTINAL ARGUMENTS
            file_name_structure    Name of file containing secondary structure and COMs information.
                                Default value: 'residues_secondaryStructure_and_COMs'
                                
            neighbor_cutoff_distance    Cutoff distance used to identify neighbor residues.
                                        If any value in the list of columns distance_columns' is beloe or equal
                                        to this distance, the residues are neighbors.
                                        Default value: 8 angstrom
            
            distance_columns    list of columns in which to search for the cutoff distances
                                Default value: None (search in all distances columns)
            
            pdb_name    structure name to use in the clened file: {pdb_name}_distance_map.csv
                        Default value: directory name
                        
            save        save cleaned file as CSV in the directory
                        Default value: True 
        """
        return clean_neighbor_residues(*arg,**kwargs)


#=====================================================
#===== END
#=====================================================
if __name__ == "__main__":
    main()