# All-electron reference data

Here we summarize the steps how to get all-electron PBE results from the previous
precision benchmark of (E. Bosoni et al., How to verify the precision of density-functional-theory implementations via reproducible and universal workflows, Nat. Rev. Phys. 6, 45 (2024), doi: 10.1038/s42254-023-00655-3).

Note that this reference data is only used in plotting routines,
and not needed to re-run the optimization itself.

## Downloading raw all-electron results

- install aiida 
  https://aiida.readthedocs.io/projects/aiida-core/en/stable/intro/install_system.html#intro-get-started-system-wide-install

- All datasets are at https://archive.materialscloud.org/record/2023.81 (can copy their link by right-click)

- If using venv: Activate venv where aiida was installed
```
source .venv/bin/activate
```
- Download e.g. unaries and oxides for Wien2k
  ```
  verdi archive import https://archive.materialscloud.org/record/file\?record_id\=1770\&filename\=acwf-verification_oxides-verification-PBE-v1_results_wien2k.aiida

  verdi archive import https://archive.materialscloud.org/record/file\?record_id\=1770\&filename\=acwf-verification_unaries-verification-PBE-v1_results_wien2k.aiida
  ```
To see the actual names of the downloaded groups, do (with capital A)
```
verdi group list -A
```
Now execute
```
python get_energy_volume_data_from_aiida.py  # Prints energy-volume curve data
```
By default this prints Li-BCC data from wien2k. For other data, adapt
the query script at get_energy_volume_data_from_aiida.py.

Only to reproduce the analysis.ipynb plotting routines, you need to write the printed
results to two separate files, `wien2k_Li-BCC.json` and `wien2k_Li-XO.json`.

