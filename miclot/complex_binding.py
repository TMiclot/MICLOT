#!/usr/bin/env python3

"""
This script is part of MICLOT it calculate binding affinity, following the method of
 the contacts-based method by [Vangone *et al.* (2015)](https://doi.org/10.7554/eLife.07454).
 This code is different and is not related to the one from
 [PRODIGY on GitHub](https://github.com/haddocking/prodigy/).
"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"

__all__ = ['identify_contacts','identify_NIS_residues_SASA','compute_binding_energy']


#=====================================================
#===== Import modules
#====================================================
import os
import mdtraj as md
import numpy as np
import pandas as pd
pd.options.mode.copy_on_write = True
import freesasa
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

# Import module from MICLOT
from miclot.utilities import pdb2pandas, mdtraj_chainID_2_chainName



#=====================================================
#===== Fucntion to calculate contact & identify contact types between residues
#====================================================
def identify_contacts(pdb_file_path, chainName_receptor, chainName_ligand, write_outfile=True):
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
#===== Function to prepare structure for freesasa
#=====================================================
class StructureFilter(Select):
    def __init__(self, allowed_chains):
        super().__init__()
        self.allowed_chains = allowed_chains

    def accept_chain(self, chain):
        # Keep only chains "C" and "G"
        return chain.get_id() in self.allowed_chains

    def accept_residue(self, residue):
        # Exclude heteroatoms (e.g., water) and keep only standard amino acids or nucleotides
        return not residue.id[0].startswith("W") and not residue.id[0].startswith("H")

    def accept_atom(self, atom):
        # Exclude hydrogens
        return atom.element != "H" and not atom.is_disordered()


#---------------------------------------------------------------------------
def filter_structure(input_pdb_file, allowed_chains=["A"], model_id=0):
    """
    Is not defined to be use directly by the user.
    Use the class StructureFilter to clean a PDB file and output the result as Biopython structure format.

    ARGUMENTS
    allowed_chains   Chain name to keep. Default is A.
    model_id         Model ID to keep. Default is 0
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('filtered_structure', input_pdb_file)

    # Select only the specified model
    model = structure[model_id]

    # Define output using PDBIO
    io = PDBIO()
    io.set_structure(model)  # Set the specific model

    # Use the custom StructureFilter to apply the filter
    filtered_structure = StructureFilter(allowed_chains)
    
    return io, filtered_structure  # Return IO object and filter for further use or saving





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
    #----- Clean structures -----
    io_complex, structure_complex   = filter_structure(pdb_file_path, chainName_complex)
    io_receptor, structure_receptor = filter_structure(pdb_file_path, chainName_receptor)
    io_ligand, structure_ligand     = filter_structure(pdb_file_path, chainName_ligand)

    #----- Save cleaned structures as PDB -----
    # Get the output path with the original PDB name
    file_path = pdb_file_path.replace('.pdb', '') # Replace '.pdb' with an empty string
    # Save PDB
    io_complex.save(f"{file_path}_selected_complex.pdb", select=structure_complex)
    io_receptor.save(f"{file_path}_selected_receptor.pdb", select=structure_receptor)
    io_ligand.save(f"{file_path}_selected_ligand.pdb", select=structure_ligand)



    #===== Compute SASA calculation with freeSASA for complex, receptor and ligand =====
    # Default parameters for freesasa:
    # print(freesasa.Parameters.defaultParameters)
    # {'algorithm': 'LeeRichards',
    #  'probe-radius': 1.4,
    #  'n-points': 100,
    #  'n-slices': 20,
    #  'n-threads': 1}
    
    # Read cleaned PDB
    parser = PDBParser(QUIET=True)
    pdb_complex  = parser.get_structure('structure_complex', f"{file_path}_selected_complex.pdb")
    pdb_receptor = parser.get_structure('structure_receptor', f"{file_path}_selected_receptor.pdb")
    pdb_ligand   = parser.get_structure('structure_ligand', f"{file_path}_selected_ligand.pdb")


    # Compute freesasa
    result_complex,  sasa_classes_complex  = freesasa.calcBioPDB(pdb_complex, classifier=freesasa.Classifier.getStandardClassifier('naccess'))
    result_receptor, sasa_classes_receptor = freesasa.calcBioPDB(pdb_receptor, classifier=freesasa.Classifier.getStandardClassifier('naccess'))
    result_ligand,   sasa_classes_ligand   = freesasa.calcBioPDB(pdb_ligand, classifier=freesasa.Classifier.getStandardClassifier('naccess'))

    # create s list to store all results from freesasa
    list_result_freesasa = [result_complex, result_receptor, result_ligand]



    #===== Generate a pandas dataframe from the PDB file and contaiing chainName, resName and resSeq =====
    # Get the PDB as Pandas dataframe
    df_pdb = pdb2pandas(pdb_file_path)
    # sometyme pandas perform strange convertion of resSeq number, this commands ensure it to the a string of int number: ex: '10'
    df_pdb['resSeq'] = df_pdb['resSeq'].astype(float)
    df_pdb['resSeq'] = df_pdb['resSeq'].astype(int)
    df_pdb['resSeq'] = df_pdb['resSeq'].astype(str)

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
                df_SASA_result.loc[(df_SASA_result[f'chainName'] == str(chain)) & (df_SASA_result[f'resSeq'] == str(resSeq)), f'{column_name}'] = residue_asa_total


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
    else:
        os.remove(f"{file_path}_selected_complex.pdb")
        os.remove(f"{file_path}_selected_receptor.pdb")
        os.remove(f"{file_path}_selected_ligand.pdb")
    

    #===== Return the final dataframe =====
    return df_SASA_result_filtered





#====================================================
#===== Function to make report on 
#====================================================
def compute_binding_energy(pdb_file_path, chainName_receptor, chainName_ligand, temperature_celsius=25, write_outfile=True):
    """
    DESCRIPTION
        Compute DeltaG and Kd, then return a final report containing all contacts and NIS informations.
    
    ARGUMENTS
        pdb_file_path         Path to the PDB file.
        chainName_receptor    List of chain names of the receptor.
        chainName_ligand      List of chain names of the ligand.

    OPTIONAL ARGUMENTS
        temperature_celsius    Temperature in Celsius.
                               Default value: 25.0

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

    # Check if the index exists
    for index_to_check in new_indices:
        if index_to_check not in series_counts_types.index:
            series_counts_types[index_to_check] = 0  # Create the index with a value of '0'

    # change indices
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

    # Check if the index exists
    for index_to_check in new_indices:
        if index_to_check not in series_SASA_types.index:
            series_SASA_types[index_to_check] = 0  # Create the index with a value of '0'

    # change indices
    series_SASA_types = series_SASA_types.rename(index=new_indices)

    # Get the sum of all residue types
    series_SASA_types['TOTAL_NIS_types'] = series_SASA_types.sum()


    # Get proportion of residue types (polar, apolr, charged) & convert it to %
    series_SASA_types_percentages = df_SASA['resType'].value_counts(normalize=True) *100

    # Renames indices
    new_indices = {'polar': 'NIS_polar(%)',
                   'apolar': 'NIS_apolar(%)',
                   'charged': 'NIS_charged(%)'}
    series_SASA_types_percentages = series_SASA_types_percentages.rename(index=new_indices)


    #----- Create the final Series for report by concatenate all the previous -----
    report = pd.concat([series_counts_types, series_SASA_types, series_SASA_types_percentages])


    #===== compute DeltaG and kd =====
    #----- Setup initial parameters -----
    temperature_kelvin  = temperature_celsius + 273.15
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
    report['temperature(C)'] = temperature_celsius
    report['temperature(K)'] = temperature_kelvin
    report['DeltaG(kcal/mol)'] = Delat_G
    report[f'dissociation_constant_(M)_at_{temperature_celsius}(C)'] = float(kd)
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
        df_report.to_csv(f'{file_path}_report_interaction.csv', index=False)
    

    #===== Return the final dataframe =====
    return df_report




#====================================================
#=====| Module end |=====
if __name__ == "__main__":
    main()
