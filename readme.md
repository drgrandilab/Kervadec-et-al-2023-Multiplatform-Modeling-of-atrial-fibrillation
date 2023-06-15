
## Multiplatform modeling of atrial fibrillation



This folder contains source code for modeling human atrial action potential and Ca2+ with changes with PLN and beta-adrenergic activation.   

---

### Single Cells:


#### Folders
* Single-Cells: everything needed for single cell simulations and analysis
* Logistic_Regression: logstic regression analysis
* EAD_DAD_Analysis: analysis code for identifying EAD/DAD from APs. 
* parameter perturbation: para_log.dat.New

#### Tools
* compiler: C++ with CXX=mpiicpc
* CVODE with Sundials 
* zlib
* population run with python (run_2Hz.py)

To run the source code: 
```bash
make clean;
make 
python run_2Hz.py
```
---

### Tissue modeling

#### Folder: 2D-Tissue-src
* lib and lib_tissue: library for tissue model
* Geometry: tissue geometry 
* Initial_Conditions: ICs for cells/nodes in tissue
* Python_Scripts: post analysis code in python

#### Tools
* compiler: C++ with CXX=mpiicpc (for MPI)
* CVODE with Sundials 
* zlib (to read in geometry data)
* population run with python (run_2Hz.py)
* slurm job management


To run the source code: 
```bash
make clean;
make 

sbatch run.slurm
# or
sh run.sh

```