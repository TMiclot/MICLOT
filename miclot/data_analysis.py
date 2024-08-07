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
    return dataframe
    
    list_neighbor = ["0_5", "3_4", "10_3", "1_9"]
    list_couple = ["0_1", "3_5", "10_12", "14_15"]

    df = find_neighbor_couples(list_neighbor, list_couple)
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




###################################################### Cleaning & Concatenate

#=====================================================
#===== Function to concatenate CSV files
#====================================================
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
        plot_minimization(*arg,**kwargs)





#=====================================================
#===== Class to concatenate (cleaned) CSV files
#=====================================================
class concatenate:
    def __init__():
        """
        Functions used to concatenate CSV files
        """
        pass
    
    #----- Concatenate CSV files over subdirectory ----- 
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
        
    #----- Concatenate files with ASA information -----
    def ASA(*arg, **kwargs):
        """
        Concatenate '*ASA_complex_receptor_ligand_and_interface.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'ASA_complex_receptor_ligand_and_interface.csv', 'final_ASA_complexes.csv', logfile_name='concatenate_ASA', **kwargs)
        
    #----- Concatenate files with interaction_table information -----
    def interaction_table(*arg, **kwargs):
        """
        Concatenate '*clean_interactions_table.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'clean_interactions_table.csv', 'final_interactions_table.csv', logfile_name='concatenate_interactions_table', **kwargs)
        
    #----- Concatenate files with neighbor_residues information -----
    def neighbor_residues(*arg, **kwargs):
        """
        Concatenate '*clean_neighbor_residues.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'clean_neighbor_residues.csv', 'final_neighbor_residues.csv', logfile_name='concatenate_neighbor_residues', **kwargs)
           
    #----- Concatenate files with neighbor_pairs information -----
    def neighbor_pairs(*arg, **kwargs):
        """
        Concatenate '*clean_neighbor_pairs.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'clean_neighbor_pairs.csv', 'final_neighbor_pairs.csv', logfile_name='concatenate_neighbor_pairs', **kwargs)
        
    #----- Concatenate files with structure information -----
    def structure(*arg, **kwargs):
        """
        Concatenate '*clean_structure.csv' files in all subdirectories in a directory.
        
        ARGUMENT
            directory    directory containing all subdirectory with CSV files.
        """
        concatenate_csv(*arg, 'clean_structure.csv', 'final_structure.csv', logfile_name='concatenate_', **kwargs)



#=====================================================
#===== END
#=====================================================
if __name__ == "__main__":
    main()