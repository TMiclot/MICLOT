#!/usr/bin/env python3

"""
This Database contain propeerties of:
    - amino acids
    - elements found in proteins
"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"


__all__ = ['amino_acids', 'elements']





#=====================================================
#===== Import modules
#=====================================================
import os
import pandas as pd





#=====================================================
#===== Fucntion to load anlys CSV file from database.py directory
#=====================================================
def load_csv_database(CSV_filename):
    # Get the directory where the database.py file is located
    script_dir = os.path.dirname(__file__)
    # Create the full path to the CSV file
    file_path = os.path.join(script_dir, CSV_filename)
    # Read and return the CSV file
    return pd.read_csv(file_path)





#=====================================================
#===== Variable containing DB as Pandas table
#=====================================================

#=====
amino_acids = load_csv_database('DB_amino_acids_properties.csv')

#=====
elements = load_csv_database('DB_elements_properties.csv')



#=====================================================
#=====| Module end |=====
if __name__ == "__main__":
    main()