#===== Initialise directories & Dataset =====
# Directory with initial pdb files
directory_dataset_structures = '/Data/analysis/FINAL_dataset_PDB'

# CSV dataset with all informations
dataset_file = '/Data/analysis/FINAL_dataset_complete.csv'

# Create working directory
directory_working = '/Data/analysis/analysis'

# Directory of MICLOT
directory_MICLOT = '/Data/analysis/MICLOT'

# Logging file
log_file_path = 'script.log'




#############################################
## python modules & Configuration

#===== Import modules =====
import sys
import os
import shutil
import glob
import logging
from time import sleep


#===== Limit number of CPU usage =====
# get process ID of this script
pid = os.getpid()

#Get total number of CPU -2 to avoid crash
max_cpu = os.cpu_count() -4

# Set the process to the range of from 0 to max_cpu
os.sched_setaffinity(pid, range(max_cpu))


#===== Import third-party modules =====
import mdtraj as md
import pandas as pd
from pdbfixer import PDBFixer
from openmm.app import PDBFile

#import matplotlib.pyplot as plt
#plt.ioff() #disables interactive mode and storing figure in a variable won't display any output



#===== Configure log file =====
logger = logging.getLogger(__name__)

# Format log output
#    2023-06-05 16:23:52,625 - INFO :: any information user must know
logging.basicConfig(format='%(asctime)s - %(levelname)s :: %(message)s',
                    filename=f'{log_file_path}',
                    encoding='utf-8',
                    level=logging.INFO)

# append log file
logger.info(">>> START of execution <<<")




#############################################
## Create a termination function
def termination(error=False):
    """
    Identify causes of script ending. Then quit python.
    """
    # Ending without errors
    if error != True:
        logger.info("> Script stop normaly")

    # Ending with errors
    elif error == True:
        logger.error("> Script stop with error")

    logger.info(">>> END of execution <<<")
    
    # Quit for jupyter notebook
    raise Exception("Jupyter END of execution.") # only for jupyter
    
    # Quit for script executed in terminal
    #sys.exit()




#############################################
## MICLOT package

#===== Check MICLOT directory =====
# Return error and stop script if the directory provided for MICLOT location don't exist
if not os.path.exists(directory_MICLOT):
    # Notify user in log file
    logger.error(f"Issue with the MICLOT directory '{directory_MICLOT}'.")
    # End script with error
    termination(True)



#===== Try to add MICLOT to sys.path =====
try:
    # Add the 'MICLOT' package directory to sys.path
    if directory_MICLOT not in sys.path:
        sys.path.insert(0, directory_MICLOT)
    
    # Notify user in log file
    logger.info("Add the MICLOT package directory to sys.path")

except:
    # Notify user in log file
    logger.error("Unable to add the MICLOT package directory to sys.path")
    # End script with error
    termination(True)



#===== Import MICLOT =====
try:
    import miclot.utilities as mcu
    import miclot.interactions as mci
    import miclot.complex_binding as mcb
    import miclot.cys_bridges as mcy
    
    # Notify user in log file
    logger.info("MICLOT is successfully imported")

except:
    # Notify user in log file
    logger.error("Unable to import MICLOT")
    # End script with error
    termination(True)



    
#############################################
## Working directory
try:
    if not os.path.exists(directory_working):
        # create the directory
        os.makedirs(directory_working)
        
        # Notify user in log file
        logger.info(f"Working directory '{directory_working}' doesn't exist, so it's been created")

    elif os.path.exists(directory_working):
        # Notify user in log file
        logger.info(f"Working directory '{directory_working}' exist, and will be used as working directory")
        
except:
    # Notify user in log file
    logger.error(f"Issue with working directory '{directory_working}'")
    # End script with error
    termination(True) 
    

    
    
#############################################    
## Read the dataset file as Pandas dataframe

try:
    df_dataset_file = pd.read_csv(dataset_file)
    # Notify user in log file
    logger.info(f"Dataset file '{dataset_file}' is loaded as dataframe")

except:
    # Notify user in log file
    logger.error(f"Issue with dataset file '{dataset_file}'")
    # End script with error
    termination(True) 

    
    
    
#############################################
## 0. Progress & Create dedicated directory & copy initial structure
## 1. Clean structure
## 2. Structure properties
## 3. Analysis

len_df_dataset_file = len(df_dataset_file)
total_step_of_analysis = int(len_df_dataset_file)# /2)

step_current = 0

# loop over all structures in the dataset
for index, row in df_dataset_file.iterrows():
    # !!! For testing use the first row
    if index > total_step_of_analysis :
        logger.info(f"Exclude structure {row['pdb']} in row index: {index}")
        continue
        
        
    #===== 0a) Notify user on progress =====
    # calculate progress
    step_current += 1
    step_percent = (step_current / total_step_of_analysis) *100
    
    # Notify user in log file
    logger.info(f"> Progress: {step_current}/{total_step_of_analysis} ({step_percent} %)")
    
    
    
    #===== 0b) Get structure to work on =====
    #----- Get PDB code -----
    pdb = row['pdb']    
    
    #----- Identify file in the working directory with the PDB code -----
    try:
        file_pattern = os.path.join(directory_dataset_structures, f'{pdb}*.pdb')
        initial_pdb_file_path = glob.glob(file_pattern)

        if 1 < len(initial_pdb_file_path):
            # Notify user in log file
            logger.warning(f"Multiple files with PDB code '{pdb}'. Skip this PDB and continue with next item. File are {initial_pdb_file_path}")
            # Go to next structure
            continue
        
        else:
            initial_pdb_file_path = initial_pdb_file_path[0]
    
    except:
        # Notify user in log file
        logger.warning(f"No file with PDB code '{pdb}'. Skip this PDB and continue with next item.")
        # Go to next structure
        continue

        
    #----- if their is structure file with PDB code -----
    # Notify user in log file
    logger.info(f"Analysis will be done on PDB code '{pdb}', with initial structure '{initial_pdb_file_path}'")
    
    
    
    
    #===== 0c) Create a dedicated directory in the working_directory =====
    #----- create path for the dedicated directory -----
    pdb_directory = f'{directory_working}/{pdb}'
    
    #----- check if the directory dont't exist -----
    if not os.path.exists(pdb_directory):
        # Create the directory
        os.makedirs(pdb_directory)
        
        # Notify user in log file
        logger.info(f"Dedicated directory is created in '{pdb_directory}'")
    
    else:
        # Notify user in log file
        logger.warning(f"Dedicated directory for PDB code '{pdb}' exist. Skip this PDB and continue with next item. Directory is {pdb_directory}")
        # Go to next structure
        continue
    
    
    
    
    #===== 0d) Copy structure in it's dedicated directory =====
    # Get proper file name
    pdb_file_name_initial = initial_pdb_file_path.split('/')[-1]
    
    # copy the file
    try:
        # copy the initial structure to a dedicated directory in the 'directory_working'
        shutil.copy(initial_pdb_file_path, pdb_directory)
        
        # Notify user in log file
        logger.info(f"Copy initial structure file '{pdb_file_name_initial}' in '{pdb_directory}'")
        
        pdb_initial = f'{pdb_directory}/{pdb_file_name_initial}'
    
    except:
        # Notify user in log file
        logger.warning(f"Unable to copy initial structure file '{pdb_file_name_initial}'. Skip this PDB and continue with next item. Directory is '{pdb_directory}'")
        # Go to next structure
        continue    
        
        
        
        
    #===== 1a) Remove HETATM lines =====
    sleep(1)
    # set name for file without hetatm
    pdb_noHETATM_name = pdb_file_name_initial.split('.')[0]
    pdb_noHETATM = f'{pdb_directory}/{pdb_noHETATM_name}_noHETATM.pdb'
    
    try:
        # loop over the file and write only line not starting with 'HETATM' in the output file
        with open(pdb_initial, 'r') as infile, open(pdb_noHETATM, 'w') as outfile:
            for line in infile:
                if not line.startswith('HETATM'):
                    outfile.write(line)
        
        # Notify user in log file
        logger.info(f"Remove HETATM in '{pdb_initial}' and write '{pdb_noHETATM}'")
    
    except:
        # Notify user in log file
        logger.warning(f"Unable to remove HETATM '{pdb_initial}'. Skip this PDB and continue with next item.")
        # Go to next structure
        continue        
                
    
    
    
    #===== 1b) Perform PDB2PQR =====
    # Add missing atoms, rename residue with their protonation state, remove water only
    sleep(1)
    try:
        # Run PDB2PQR
        mcu.pbd2pqr_parse(pdb_noHETATM, ph=row['pH']) # use pH value as define in the dataset
        
        # get generated file path
        file_pattern = os.path.join(pdb_directory, f'*_pqr.pdb')
        pdb_PDB2PQR = glob.glob(file_pattern)[0]

        # Notify user in log file
        logger.info(f"Perform PDB2PQR on '{pdb_noHETATM}' and export '{pdb_PDB2PQR}'")
    
    except:
        # Notify user in log file
        logger.warning(f"Unable to perform PDB2PQR on '{pdb_noHETATM}'. Skip this PDB and continue with next item.")
        # Go to next structure
        continue
                
    
    
    
    #===== 1c) Fix output of PDB2PQR from  =====
    ## ! This part rename residues with their standard names !
    sleep(1)
    try:
        # use PDBfixer to fix pdb file for openmm
        pdbfixer = PDBFixer(filename=pdb_PDB2PQR)
        pdbfixer.removeHeterogens(True) # Ensure remove heterogens
        ##pdbfixer.findNonstandardResidues()
        ##pdbfixer.replaceNonstandardResidues() # replace nonstandard residues with their standard versions
        pdbfixer.findMissingResidues()
        pdbfixer.missingResidues = {}
        pdbfixer.findMissingAtoms()
        pdbfixer.addMissingAtoms()
        pdbfixer.addMissingHydrogens(row['pH']) # The argument is the pH based on which to select protonation states.

        pdb_FIXED_name = pdb_PDB2PQR.split('.')[0]
        pdb_FIXED = f'{pdb_FIXED_name}_fixed.pdb'

        # write fixed PDB file
        PDBFile.writeFile(pdbfixer.topology, pdbfixer.positions, open(f'{pdb_FIXED}', 'w'))

        # Notify user in log file
        logger.info(f"Fix structure file for openMM. Structure '{pdb_PDB2PQR}'")
        
    except:
        # Notify user in log file
        logger.warning(f"Unable to fix {pdb_PDB2PQR}. Skip this PDB and continue with next item.")
        # Go to next structure
        continue

    
    
    
    #===== 1d) Perform vaccum minimization =====
    sleep(1)
    try:
        # max_iterations=0 to minimize the structure until the results converge regardless the number of iterations it takes.
        mcu.minimize_pdb(pdb_FIXED, max_iterations=0)

        # get generated file path
        file_pattern = os.path.join(pdb_directory, f'*_minimized.pdb')
        pdb_MINIMIZE = glob.glob(file_pattern)[0]

        # Notify user in log file
        logger.info(f"Perform minimization on '{pdb_FIXED}' and export '{pdb_MINIMIZE}'")
        
    except:
        # Notify user in log file
        logger.warning(f"Unable to minimize {pdb_FIXED}. Skip this PDB and continue with next item.")
        # Go to next structure
        continue       
    
    
    
    #===== 1f) Fix topology bonds =====
    sleep(1)
    try:
        # fix topology with missing bonds
        mcu.fix_topology_bonds(pdb_MINIMIZE)

        # get generated file path
        file_pattern = os.path.join(pdb_directory, f'*_corrected_bonds.pdb')
        pdb_fixed_bonds = glob.glob(file_pattern)[0]

        # Notify user in log file
        logger.info(f"Fix missing bonds in the '{pdb_MINIMIZE}' and save it as {pdb_fixed_bonds}")
        
    except:
        # Notify user in log file
        logger.warning(f"Unable to fix bonds in '{pdb_MINIMIZE}'. Skip this PDB and continue with next item.")
        # Go to next structure
        continue    
    
       
    
    
    
    #===== 1g) Check for Cys-Cys bridges =====
    sleep(1)
    try:
        # Identify CYS bridges types
        mcy.cys_bridges_inPDB(pdb_fixed_bonds)
        
        # get generated file path
        file_pattern = os.path.join(pdb_directory, f'*_CYS_bridges.pdb')
        pdb_CYSbridges = glob.glob(file_pattern)[0]
        
        # Notify user in log file
        logger.info(f"Identify disulfide, diselenide or selenosulfide bridges on '{pdb_fixed_bonds}' and export '{pdb_CYSbridges}'")
        #
    except:
        # Notify user in log file
        logger.warning(f"Unable to identify disulfide, diselenide or selenosulfide bridges {pdb_fixed_bonds}. Skip this PDB and continue with next item.")
        
        
        
        
        
    #===== 2a) Get secondary structure =====
    sleep(1)
    try:
        # identify secondary structure
        mcu.get_sequence_secstruct(pdb_CYSbridges)

        # Notify user in log file
        logger.info(f"Identify secondary structure in '{pdb_CYSbridges}'")
        #
    except:
        # Notify user in log file
        logger.warning(f"Unable to identify secondary structure in '{pdb_CYSbridges}. Skip this PDB and continue with next item.")
        # Go to next structure
        continue
        
        
        
        
    #===== 2b) Get list of interacting chain as MDTraj IDs =====
    # no need to sleep 'pdb_PDB2PQR' already exist for a ""long time""
    try:
        #----- Get chains as MDTraj IDs from the 'pdb_PDB2PQR' -----
        # it use 'pdb_PDB2PQR' because minimization rename all chains, so their no correlation can be done
        # between reported interacting chains and ID
        dict_chainID_2_chainName, dict_chainName_2_chainID = mcu.mdtraj_chainID_2_chainName(pdb_PDB2PQR)

        # Notify user in log file 
        logger.info(f"Correlate chain name in '{pdb_PDB2PQR}' with corresponding MDTraj IDs")


        #----- Get interacting chain IDs -----
        # Read reported interacting chains from dataset
        interacting_chain_1 = row['chain_1']
        interacting_chain_2 = row['chain_2']

        # Identify if interacting chains are alphabet or digit (digit only for dimer comming from PDBind)
        if interacting_chain_1.isdigit() and interacting_chain_2.isdigit():

            # convert to list with integer
            list_chain_1_ID = [int(i) for i in interacting_chain_1]
            list_chain_2_ID = [int(i) for i in interacting_chain_2]

            # Notify user in log file 
            logger.info(f"Get interacting chains as MDTraj IDs for {pdb}")


        elif interacting_chain_1.isalpha() and interacting_chain_2.isalpha():

            # create list with the list_chain_1 = list(interacting_chain_1)
            list_chain_1 = list(interacting_chain_1)
            list_chain_2 = list(interacting_chain_2)

            # Create the list of chain IDs
            list_chain_1_ID = sum([dict_chainName_2_chainID[i] for i in list_chain_1], [])
            list_chain_2_ID = sum([dict_chainName_2_chainID[i] for i in list_chain_2], [])

            # Notify user in log file 
            logger.info(f"Get interacting chains as MDTraj IDs for {pdb}")

        else:
            # Notify user in log file
            logger.warning(f"Unable to differentiate between alphabet and digit for interacting chains in the dataset for PDB {pdb}. Skip this PDB and continue with next item.")
            # Go to next structure
            continue
        
  
    except:
        # Notify user in log file
        logger.warning(f"Unable to get chains as MDTraj IDs for '{pdb}' from dataset. Skip this PDB and continue with next item.")
        # Go to next structure
        continue
    
    
    
    
    #===== 2c) Identify protein region of residues =====
    # like previously no need to sleep
    try:
        # identify protein region of all residues
        mcu.get_protein_region(pdb_CYSbridges, list_chain_1_ID, list_chain_2_ID)

        # Notify user in log file
        logger.info(f"Identify protein region for all residues in '{pdb_CYSbridges}'")
        
    except:
        # Notify user in log file
        logger.warning(f"Unable to identify protein region for all residues in '{pdb_CYSbridges}'. Skip this PDB and continue with next item.")
        # Go to next structure
        continue
        
        
        
        
    #===== 3) Identify all interaction types between residues =====
    sleep(1)
    try:
        # Get table file name
        table_file_name = pdb_CYSbridges.split('.')[0]

        # Notify user in log file
        logger.info(f"Make table of interaction for '{pdb_CYSbridges}' -- Start ")

        # Load molecular structure using MDTraj
        traj = md.load(pdb_CYSbridges, top=pdb_CYSbridges)

        # identify protein region of all residues
        mci.interaction_table_whole_system(traj, path_table_outfile=f"{table_file_name}_interaction_table_whole_system.csv", path_class_outfile=f"{table_file_name}_class")

        # Notify user in log file
        logger.info(f"Make table of interaction for '{pdb_CYSbridges}' -- End ")
        
    except:
        # Notify user in log file
        logger.warning(f"Unable to make table of interaction for '{pdb_CYSbridges}'. Skip this PDB and continue with next item.")
        # Go to next structure
        continue

    
    
#############################################
## End of the script
termination()
