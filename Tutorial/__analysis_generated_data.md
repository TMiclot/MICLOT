[Tutorials home](Tutorials.md)

# Analyze generated data

## First step: load required modules

```python
# Define the path to of the package
package_dir = '/path/to/MICLOT'

# Add the 'miclot' package directory to sys.path
if package_dir not in sys.path:
    sys.path.insert(0, package_dir)
```

Now you can import the packages:

```python
import miclot.analysis as mca
```

## 1. Plot minimization energies

### 1.1. Normal usage

It is possible to prive directory or file path in the command.

```python
# directory path
a = plot_minimization('/Data/1acb/')
a.show()
```

```python
# file path
a = plot_minimization('/Data/1acb/1acb_noHETATM_pqr_fixed_minimization_log.csv')
a.show()
```

### 1.2. Change the structure name in the title

If you want to change the name in the title, you can modify it using `pdb_name` argument.

```python
a = plot_minimization('/Data/1acb/1acb_noHETATM_pqr_fixed_minimization_log.csv', pdb_name='Hydrolase/Hydrolase inhibitor')
a.show()
```

### 1.3. Saving the graph

```python
plot_minimization('/Data/1acb/1acb_noHETATM_pqr_fixed_minimization_log.csv', save_graph=True)
```

### 1.4. Example of output

<img src="tuto_pictures/minimization_energies.png" width="900">
