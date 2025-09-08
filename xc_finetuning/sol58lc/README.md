- The original 58 solids were collected in https://journals.aps.org/prb/abstract/10.1103/PhysRevB.79.085104
- The 58 experimental lattice constants **including zero-point corrections** (based on PBE & DFPT in Quantum Espresso), are from Table III https://journals.aps.org/prb/abstract/10.1103/PhysRevB.85.014111
- the CSV file itself is from the supplement of mBEEF-vdW 2016 https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.235162


Here we use the atomic simulation environment (ASE) to generate `.extxyz` files.


The files can be read in julia as follows:
```
using AtomsIO
system = load_system("GaAs_b3.extxyz")
```
