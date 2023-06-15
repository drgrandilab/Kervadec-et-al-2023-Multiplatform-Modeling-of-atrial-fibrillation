
## Multiplatform modeling of atrial fibrillation



This folder contains source code for modeling human atrial action potential and Ca2+ with changes with PLN and beta-adrenergic activation.   

Please cite the following publication when applicable:

* Kervadec A, Kezos J, Ni H, Yu M, Marchant J, Spiering S, Kannan S, Kwon C, Andersen P, Bodmer R, Grandi E, Ocorr K, Colas AR. Multiplatform modeling of atrial fibrillation identifies phospholamban as central regulator of cardiac rhythm. Dis Model Mech. 2023 Jun 9;. doi: 10.1242/dmm.049962. [Epub ahead of print] PubMed PMID: 37293707. 
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
