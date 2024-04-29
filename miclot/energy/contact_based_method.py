#!/usr/bin/env python3

"""
This script is part of MICLOT ...
"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"


#=====================================================
#===== Import modules
#====================================================
import mdtraj as md
import numpy as np
import pandas as pd
pd.options.mode.copy_on_write = True
import freesasa
from Bio.PDB import PDBParser, Polypeptide

# Import module from MICLOT
from ../utilities import pdb2pandas, mdtraj_chainID_2_chainName






#=====================================================
#===== Fucntion to calculate contact & identify contact types between residues
#====================================================
def identify_contacts(pdb_file_path, chainName_receptor, chainName_ligand, write_logfile=True):
    """
    DESCRIPTION
        Use MDTraj to calculate contacts types between residues, based on their closest heavy atom distance.
        The cutoff used is 5.5 angstroms.
        Contacts types are:
            - apolar-polar
            - apolar-apolar
            - charged-polar
            - polar-polar
            - apolar-charged
            - charged-charged
            
    RETURN
        pandas dataframe with residues in contacts, with:
            chainIDs, resIDs, resNames, distance, contact_type
            
    ARGUMENTS
        pdb_file_path       Path to the pdb file (string format).
        chainName_receptor    chain ID of the receceptor list of string format.
                            Example: ['A','B']
        chainName_ligand      chain ID of the ligand in list of string format.
                            Example: ['C']

    OPTIONAL ARGUMENTS
        write_outfile       (True/False) write the output CSV files.
                            Default vaule: True
    """
    
    #===== Dictionnary to store AA properties when consider number of contact =====
    dict_amino_acid_properties_contacts = {
        "ALA": "apolar",
        "CYS": "apolar",
        "GLU": "charged",
        "ASP": "charged",
        "GLY": "apolar",
        "PHE": "apolar",
        "ILE": "apolar",
        "HIS": "charged",
        "LYS": "charged",
        "MET": "apolar",
        "LEU": "apolar",
        "ASN": "polar",
        "GLN": "polar",
        "PRO": "apolar",
        "SER": "polar",
        "ARG": "charged",
        "THR": "polar",
        "TRP": "apolar",
        "VAL": "apolar",
        "TYR": "apolar",
        }
    
    #===== Load the PDB file with MDTraj=====
    traj = md.load(pdb_file_path, top=pdb_file_path)

    # Get dictionnaties corresponding MDTraj chainID of the receptor and ligand
    chainID_2_chainName, chainName_2_chainID = mdtraj_chainID_2_chainName(pdb_file_path, write_outfile=False)

    # Using list comprehension to populate with the dictionnarie 'chainName_2_chainID' to convert all chainName into their chainID
    #     Note that in MDTraj one chainName can have multiple chainIDs
    chainID_receptor = []
    chainID_receptor.extend([chainID for chain in chainName_receptor for chainID in chainName_2_chainID[chain]])
    chainID_ligand = []
    chainID_ligand.extend([chainID for chain in chainName_ligand for chainID in chainName_2_chainID[chain]])

    #===== Get contact between residues in receptor and ligand =====
    # Create list to store results
    list_residues_ID_receptor = []
    list_residues_ID_ligand = []
    list_residues_name_receptor = []
    list_residues_name_ligand = []
    list_residues_chainID_receptor = []
    list_residues_chainID_ligand = []
    list_contact_distance = []
    list_contact_type = []



    #===== Get all CA atom index in receptor or ligand =====
    # ensure chainID are string
    str_chainID_receptor = [str(i) for i in chainID_receptor]
    str_chainID_ligand = [str(i) for i in chainID_ligand]

    # join list of chainID as sigle string
    str_chainID_receptor = ' '.join(str_chainID_receptor)
    str_chainID_ligand = ' '.join(str_chainID_ligand)

    # Get the ID of all CA atoms in receptor or ligand
    select_CA_receptor = traj.topology.select(f"chainid {str_chainID_receptor} and protein and name CA")
    select_CA_ligand = traj.topology.select(f"chainid {str_chainID_ligand} and protein and name CA")



    #===== Search all contact between residue in receptor and ligand =====
    # loop over all CA atoms in receptor and ligand
    for CA_receptor in select_CA_receptor:
        for CA_ligand in select_CA_ligand:

            # Get the residue of receptor (_R), or ligand (_L), corresponding to CA atom
            residue_R = traj.topology.atom(CA_receptor).residue
            residue_L = traj.topology.atom(CA_ligand).residue

            # Get distance between the two residues
            distance = float(md.compute_contacts(traj, [[residue_R.index, residue_L.index]], scheme='closest-heavy')[0][0] *10)

            # select only pair with contact distance <=5.5 angstrom
            if (distance <= 5.5) and (residue_R.name in dict_amino_acid_properties_contacts) and (residue_L.name in dict_amino_acid_properties_contacts):

                # Append list of distance
                list_contact_distance.append(distance)

                # Append list of contact type
                type_residue_R = dict_amino_acid_properties_contacts[residue_R.name]
                type_residue_L = dict_amino_acid_properties_contacts[residue_L.name]
                contact_type = '-'.join( sorted([type_residue_R, type_residue_L]) )
                list_contact_type.append(contact_type)

                # Get residues Receptor properties and append lists
                list_residues_ID_receptor.append(residue_R.index)
                list_residues_name_receptor.append(residue_R.name)
                list_residues_chainID_receptor.append(residue_R.chain.index)

                # Get residues Ligand properties and append lists
                list_residues_ID_ligand.append(residue_L.index)
                list_residues_name_ligand.append(residue_L.name)
                list_residues_chainID_ligand.append(residue_L.chain.index)


    # Create a pandas dataframe with the result of contacts            
    df_contact = pd.DataFrame({'residues_chainID_receptor':list_residues_chainID_receptor,
                              'residues_ID_receptor':list_residues_ID_receptor,
                              'residues_name_receptor':list_residues_name_receptor,
                              'residues_chainID_ligand':list_residues_chainID_ligand,
                              'residues_ID_ligand':list_residues_ID_ligand,
                              'residues_name_ligand':list_residues_name_ligand,        
                              'contact_distance':list_contact_distance,
                              'contact_type':list_contact_type,
                              })


    #===== Save output files =====
    if write_outfile == True:
        # Get the output path with the original PDB name
        file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string

        # save the concatenate dataframe to a file
        df_contact.to_csv(f'{file_path}_residue_contact_types.csv', index=False)


    #===== return final dataframe =====
    return df_contact






#=====================================================
#===== Function to identify NIS(interface) residues based on their SASA
#====================================================
def identify_NIS_residues_SASA(pdb_file_path, chainName_receptor, chainName_ligand, write_outfile=True):
    """
    DESCRIPTION
        Return all residues corresponding to NIS/interface

    ARGUMENTS
        pdb_file_path       Path to the pdb file (string format).
        chainName_receptor    chain ID of the receceptor list of string format.
                            Example: ['A','B']
        chainName_ligand      chain ID of the ligand in list of string format.
                            Example: ['C']

    OPTIONAL ARGUMENTS
        write_outfile       (True/False) write the output CSV files.
                            Default vaule: True   
    """
    #===== Create dictionnaries to stores AA properties =====
    # Max ASA values from NACCESS software (http://www.bioinf.manchester.ac.uk/naccess/)
    NACCESS_maxASA = {"ALA": 107.95,
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

    # Dictionnary to store AA properties when consider interface
    # it have diffrences with dict_amino_acid_properties_contacts
    dict_amino_acid_properties_interface = {
        "ALA": "apolar",
        "CYS": "polar",
        "GLU": "charged",
        "ASP": "charged",
        "GLY": "apolar",
        "PHE": "apolar",
        "ILE": "apolar",
        "HIS": "polar",
        "LYS": "charged",
        "MET": "apolar",
        "LEU": "apolar",
        "ASN": "polar",
        "GLN": "polar",
        "PRO": "apolar",
        "SER": "polar",
        "ARG": "charged",
        "THR": "polar",
        "TRP": "polar",
        "VAL": "apolar",
        "TYR": "polar",
    }
    
    #===== Setup complex chains =====
    chainName_complex = chainName_receptor + chainName_ligand


    #===== Use biopython to prepare strcture for freeSASA =====
    structure_complex   = freesasa_structure_preparation(pdb_file_path, chainName_complex)
    structure_receptor  = freesasa_structure_preparation(pdb_file_path, chainName_receptor)
    structure_ligand    = freesasa_structure_preparation(pdb_file_path, chainName_ligand)



    #===== Compute SASA calculation with freeSASA for complex, receptor and ligand =====
    # Default parameters for freesasa:
    # print(freesasa.Parameters.defaultParameters)
    # {'algorithm': 'LeeRichards',
    #  'probe-radius': 1.4,
    #  'n-points': 100,
    #  'n-slices': 20,
    #  'n-threads': 1}
    
    
    # Compute freesasa
    result_complex,  sasa_classes_complex  = freesasa.calcBioPDB(structure_complex, classifier=freesasa.Classifier.getStandardClassifier('naccess'))
    result_receptor, sasa_classes_receptor = freesasa.calcBioPDB(structure_receptor, classifier=freesasa.Classifier.getStandardClassifier('naccess'))
    result_ligand,   sasa_classes_ligand   = freesasa.calcBioPDB(structure_ligand, classifier=freesasa.Classifier.getStandardClassifier('naccess'))

    # create s list to store all results from freesasa
    list_result_freesasa = [result_complex, result_receptor, result_ligand]



    #===== Generate a pandas dataframe from the PDB file and contaiing chainName, resName and resSeq =====
    # Get the PDB as Pandas dataframe
    df_pdb = pdb2pandas(pdb_file_path)

    # Get only raw with CA atom and residue in the dict_amino_acid_properties_interface
    df_pdb_filtered = df_pdb[(df_pdb['name'].isin(['CA'])) & (df_pdb['resName'].isin(dict_amino_acid_properties_interface))]

    # Select only column with resName, resSeq, and chainName
    df_SASA_result = df_pdb_filtered[['chainName', 'resSeq', 'resName']]

    # Add resType to the dataframe
    df_SASA_result['resType'] = df_SASA_result['resName'].map(dict_amino_acid_properties_interface)#dict_amino_acid_properties_interface || dict_amino_acid_properties_contacts

    # Prepare SASA_bound and SASA_free columns
    df_SASA_result['SASA_bound'] = np.NaN
    df_SASA_result['SASA_free'] = np.NaN



    #===== Write SASA of resides when bound and unbound into the dataframe =====
    for result in list_result_freesasa:
        # Get the name of the column to modify 
        if result == result_complex:
            column_name = 'SASA_bound'
        else:
            column_name = 'SASA_free'

        # Get SASA for all residues including relative areas if available for the classifier used. 
        residueAreas = result.residueAreas()

        # Loop over all chain and residue in the residueAreas, then get the total ASA of the residue
        for chain in residueAreas:
            for resSeq in residueAreas[chain]:

                # Get all ASA informations of the residue
                residue_asa = residueAreas[chain][resSeq]
                residue_asa_total = residue_asa.total

                # Modify SASA_x column with the ASA value for the correspondind residue
                df_SASA_result.loc[(df_SASA_result[f'chainName'] == chain) & (df_SASA_result[f'resSeq'] == resSeq), f'{column_name}'] = residue_asa_total



    #===== Calculate rASA and drASA =====
    # calculate rSASA in complex or free using the MaxASA value in NACCESS_maxASA
    df_SASA_result['rASA_complex'] = df_SASA_result[f'SASA_bound'] / df_SASA_result['resName'].map(NACCESS_maxASA)
    df_SASA_result['rASA_free']    = df_SASA_result[f'SASA_free']  / df_SASA_result['resName'].map(NACCESS_maxASA)

    # calculate drASA = rASA_free - rASA_complex
    df_SASA_result['drASA'] = df_SASA_result[f'rASA_free'] - df_SASA_result[f'rASA_complex']


    #===== Select interface residues and return the result =====
    # Based on the code of 'analyse_nis' in https://github.com/haddocking/prodigy/blob/main/src/prodigy_prot/predict_IC.py filter is done on 'rASA_complex' (and not 'drASA')
    # and take value equal or greater than 5% (named NIS but corresponding to the interface)
    df_SASA_result_filtered = df_SASA_result[(5/100.0 <= df_SASA_result['rASA_complex'])]
   
   
    #===== Save output files =====
    if write_outfile == True:
        # Get the output path with the original PDB name
        file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string

        # save the concatenate dataframe to a file
        df_SASA_result_filtered.to_csv(f'{file_path}_residue_NIS_types.csv', index=False)
    

    #===== Return the final dataframe =====
    return df_SASA_result_filtered





#=====================================================
#===== Function to make report on 
#====================================================
def compute_energy(pdb_file_path, chainName_receptor, chainName_ligand, temperature_celcius=25, write_outfile=True):
    """
    DESCRIPTION
        Compute DeltaG and kd, then return a final report containing all contacts and NIS informations.
    
    ARGUMENTS
    """
    #===== Identify residues and their contacts =====
    df_contact = identify_contacts(pdb_file_path, chainName_receptor, chainName_ligand, write_outfile=write_outfile)

    #===== Identify NIS/interface residues =====
    df_SASA = identify_NIS_residues_SASA(pdb_file_path, chainName_receptor, chainName_ligand, write_outfile=write_outfile)


    #===== Make final series report =====
    #----- Create a serie for contact types -----
    # Count all 'contact_type'in the dataframe
    series_counts_types = df_contact['contact_type'].value_counts(normalize=False)

    # Renames indeces
    new_indices = {'apolar-polar':'contacts_apolar-polar',
                'apolar-apolar':'contacts_apolar-apolar',
                'charged-polar':'contacts_charged-polar',
                'polar-polar':'contacts_polar-polar',
                'apolar-charged':'contacts_apolar-charged',
                'charged-charged':'contacts_charged-charged',
                }
    series_counts_types = series_counts_types.rename(index=new_indices)

    # Get the sum of all contact
    series_counts_types['TOTAL_contacts'] = series_counts_types.sum()


    #----- Create serie for SASA (interface-NIS) -----
    # Get count residue types (polar, apolr, charged)
    series_SASA_types = df_SASA['resType'].value_counts(normalize=False)

    # Renames indeces
    new_indices = {'polar': 'NIS_polar',
                'apolar': 'NIS_apolar',
                'charged': 'NIS_charged'}
    series_SASA_types = series_SASA_types.rename(index=new_indices)

    # Get the sum of all residue types
    series_SASA_types['TOTAL_NIS_types'] = series_SASA_types.sum()

    # Get proportion of residue types (polar, apolr, charged) & convert it to %
    series_SASA_types_percentages = df_SASA['resType'].value_counts(normalize=True) *100

    # Renames indeces
    new_indices = {'polar': 'NIS_polar(%)',
                'apolar': 'NIS_apolar(%)',
                'charged': 'NIS_charged(%)'}
    series_SASA_types_percentages = series_SASA_types_percentages.rename(index=new_indices)


    #----- Create the final Series for report by concatenate all the previous -----
    report = pd.concat([series_counts_types, series_SASA_types, series_SASA_types_percentages])



    #===== compute DeltaG and kd =====
    #----- Setup initial parameters -----
    temperature_kelvin  = temperature_celcius + 273.15
    R = 0.0019858775 # ideal gas constant in kcal/mol

    #----- Calculate GeltaG using IC-NIS model -----
    Delat_G = ( - 0.09459 * report['contacts_charged-charged']
                - 0.10007 * report['contacts_apolar-charged']
                + 0.19577 * report['contacts_polar-polar']
                - 0.22671 * report['contacts_apolar-polar']
                + 0.18681 * report['NIS_apolar(%)']
                + 0.13810 * report['NIS_charged(%)']
                - 15.9433 ) # end of the equation
        
    #----- Calculate kd fron the DeltaG -----
    kd = np.exp(Delat_G / (R * temperature_kelvin) )


    #===== Add final informations on the report =====
    report['temperature(C)'] = temperature_celcius
    report['temperature(K)'] = temperature_kelvin
    report['DeltaG(kcal/mol)'] = Delat_G
    report[f'dissociation_constant_(M)_at_{temperature_celcius}(C)'] = float(kd)
    report['PDB_file'] = pdb_file_path
    report['chains_receptor'] = ' '.join(chainName_receptor)
    report['chains_ligand'] = ' '.join(chainName_ligand)

    # Convert the report serie into Dataframe
    df_report = report.reset_index()
    df_report.columns = ['name', 'value']


    #===== Save output files =====
    if write_outfile == True:
        # Get the output path with the original PDB name
        file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string

        # save the concatenate dataframe to a file
        df_report.to_csv(f'{file_path}_residue_NIS_types.csv', index=False)
    

    #===== Return the final dataframe =====
    return df_report