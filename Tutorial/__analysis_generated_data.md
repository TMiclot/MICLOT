[Tutorials home](Tutorials.md)

# Analyze generated data

## First step: load required modules

```python
import sys

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
a = plot.minimization('/Data/1acb/')
a.show()
```

```python
# file path
a = plot.minimization('/Data/1acb/1acb_noHETATM_pqr_fixed_minimization_log.csv')
a.show()
```

### 1.2. Change the structure name in the title

If you want to change the name in the title, you can modify it using `pdb_name` argument.

```python
a = plot.minimization('/Data/1acb/1acb_noHETATM_pqr_fixed_minimization_log.csv', pdb_name='Hydrolase/Hydrolase inhibitor')
a.show()
```

### 1.3. Saving the graph

```python
plot.minimization('/Data/1acb/1acb_noHETATM_pqr_fixed_minimization_log.csv', save_graph=True)
```

### 1.4. Example of output

<img src="tuto_pictures/minimization_energies.png" width="1000">





## C. Remove files

### C.1. Remove all generated plots

**Code**

```python
removed_files, error = mca.remove('/Data/1acb').plots
print(removed_files)
print(error)
```

**Result**

```
['/Data/1acb/1acb_noHETATM_pqr_fixed_minimization_log_plot_minimization_energies.png']
[]
```


### C.2. Remove all generated clean data

**Code**

```python
mca.remove.clean_data('/Data/1acb')
```

**Result**

Remove all CSV files with 'clean' in their name.



### C.3. Remove all generated saved class files

**Code**

```python
mca.remove.pickles('/Data/1acb')
```

**Result**

Remove all pkl.gz files with 'class' in their name.




### C.4. Remove any file(s)

**Code**

```python
mca.remove.files('/Data/1acb', name='hydrolase', file_format='pdb')
```

**Result**

Remove any pdb file containing 'hydrolase' in their name in the directory '/Data/1acb'.