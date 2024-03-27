[User Guide home](Manual.md)
# Properties of elements found in amino acids

| Element  | Symbol | Atomic weight | van der Waals radii (Å) | Covalent radii (Å) | Spin | Mag. moment | I Ground state ionization energy (eV) | II Ground state ionization energy (eV) |
| -------- | ------ | ------------- | ------------------- | -------------- | ---- | ----------- | -------------------------------- | --------------------------------- |
| Carbon   | C      | 12.011        | 1.70                | 0.75           | 0    | None        | 11.26030                         | 24.3833  |
| Hydrogen | H      | 1.008         | 1.10                | 0.32           | 1/2  | +2.79284    | 13.598433                        |   None   |
| Nitrogen | N      | 14.007        | 1.55                | 0.71           | 1    | +0.40376    | 14.5341                          | 29.6013  |
| Oxygen   | O      | 15.999        | 1.52                | 0.64           | 0    | None        | 13.61805                         | 35.1211  |
| Selenium | Se     | 78.971        | 1.90                | 1.18           | 0    | None        | 9.75239                          | 21.19    |
| Sulfure  | S      | 32.076        | 1.80                | 1.04           | 0    | None        | 10.36001                         | 23.33788 |

> [!IMPORTANT]
> Concerning software implementation it is necessary to remark:
> - In [MDTraj shrake_rupley](https://github.com/mdtraj/mdtraj/blob/master/mdtraj/geometry/sasa.py) the van der Waals raddi is diffrent for the hydrogen and is set at 1.20 Å.
> - Compare to MDTraj, the software NACCESS use van der Waals radii from [Chothia *et al.* (1976)](https://doi.org/10.1016/0022-2836(76)90191-1): O is 1.40 Å, N is 1.65 Å (trigonal) or 1.50 Å (tetrahedral), C is 1.87 Å (tetrahedral) or 1.76 (trigonal), S is 1.85 Å. These values are greater because they *include* van der Waals radii of hydrogen into non-hydrogen atoms, and selenium is missing.


### References
- Mantina, M., Valero, R., Cramer, C. J. & Truhlar, D. G. Atomic radii of the elements. in *CRC Handbook of chemistry and physics* (eds. Haynes, W. M., Lide, D. R. & Bruno, T. J.) 9–57 (CRC Press, 2016). [https://doi.org/10.1201/9781315380476](https://doi.org/10.1201/9781315380476)
- Prohaska, T. et al. Standard atomic weights of the elements 2021 (IUPAC Technical Report). *Pure and Applied Chemistry* 94, 573–600 (2022). [https://doi.org/10.1515/pac-2019-0603](https://doi.org/10.1515/pac-2019-0603)
- NIST: Basic Atomic Spectroscopic Data. [https://dx.doi.org/10.18434/T4FW23](https://dx.doi.org/10.18434/T4FW23)
- Chothia, C. The nature of the accessible and buried surfaces in proteins. *Journal of Molecular Biology* 105, 1–12 (1976). [https://doi.org/10.1016/0022-2836(76)90191-1](https://doi.org/10.1016/0022-2836(76)90191-1)
