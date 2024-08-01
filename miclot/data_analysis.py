# functions to generate clened data, final, tables, and analysis graphs, etc
# functions to delete generated csv (of cleaned data) and png files




#    #===== 1e) Plot system energy =====
#    sleep(1)
#    # It is not essential, so if cannot be not, the script continue analysis
#    try:
#        #Get energy log file
#        file_pattern = os.path.join(pdb_directory, f'*_minimization_log.csv')
#        minimization_log = glob.glob(file_pattern)[0]
#        
#        # Sample DataFrame creation for illustration purposes
#        df = pd.read_csv(f'{minimization_log}')
#
#        # Plot the "system_energy" column
#        plt.figure(figsize=(10, 6))
#        plt.plot(df['system_energy'], linewidth=2.5)
#
#        # Adding titles and labels
#        plot_working_file = pdb_PDB2PQR.split('/')[-1]
#        plt.title(f'System energy during minimization - ({plot_working_file})')
#        plt.xlabel('Iteration')
#        plt.ylabel('Energy (kJ/mol)')
#
#        # Adding a grid and legend
#        plt.grid(True)
#
#        # Save the plot as a high-quality PNG file
#        plot_name = minimization_log.split('_log.')[0]
#        plt.savefig(f'{plot_name}_system_energy.png', dpi=300)
#    
#        # Notify user in log file
#        logger.info(f"Plot of minimization system energy is saved as '{plot_name}_system_energy.png'")
#    
#    except:
#        # Notify user in log file
#        logger.debug(f"Unable to plot minimization system energy. Skip this step and continue analysis.")
 