# functions to generate clened data, final, tables, and analysis graphs, etc
# functions to delete generated csv (of cleaned data) and png files
# transformer le 'script_analysis.py' en une fonction de MICLOT ?



#=====================================================
#===== Function to plot minimization energies
#=====================================================
def plot_minimization(path, pdb_name=None, fig_size=[15, 4], linewidth=2, save_graph=False, dpi=300):
    """
    DESCRIPTION
        Function used to plot minimisation energies from the CSV output file.

    ARGUMENTS
        path      directory where to find the '*_minimization_log.csv' file.
                  Alternatively it can be also the path of the file.

    OPTIONAL ARGUMENTS
        pdb_name    By defaut the directory name where the CSV file is located.
                    Example 'my super PDB'
                    Default value: None
                  
        fig_size    list contiaing figure size [x_size, y_sizer]
                    Defaulůt value: [15,4]
        
        linewidth    line width of the graph.
                     Default value: 2
                   
        save_graph    save the graph in the same directory as the CSV file. (True/False)
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
    if pdb_name == None:
        if path.endswith(".csv"): # if file path is provided with file name
            pdb_name = path.split('/')[-2]
        else: # if directory path is provided (no file name)
            pdb_name = path.split('/')[-1]
        
    
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
    fig.suptitle(f'Energy variation during minimization - {pdb_name}', fontsize=20)

    #===== Save the plot as a high-quality PNG file =====
    if save_graph == True:
        plot_name = minimization_log.split('.')[0]
        plt.savefig(f'{plot_name}_minimization_energies.png', dpi=dpi)
    
    #===== return the plot =====
    return plt