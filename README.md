# PYF
The scripts for biosynthesis-growth relationship simulation, production simulation and metabolic engineering strategy optimization are in a Pyomo modelling framework.  

### About
The pipeline was written and tested on Windows Jupyter notebook. The core libraries essential for the pipeline including: cobra, Pyomo, Gurobi and related packages.  

### Environment
The scripts were written and tested with Python 3.6.5  
The core libraries essential for the pipeline including: Cobrapy toolkit: version --0.13.3, Pyomo package: version --6.4.0, Gurobi solver: version --10.0.1  

### Installation
The packages used to run the pipeline was listed in requirements.txt.  
To create a stand-alone environment named PYF with Python 3.6.5 and all the required package versions (especially for cobrapy is also available), run the following code:

1. create PYF environment using conda:  
$ conda create -n PYF python=3.6.5  

2. install related packages using pip:  
$ conda activate PYF  
$ pip install cobra==0.13.3  
$ pip install -r requirements.txt  
$ python -m ipykernel install --user --name PYF --display-name "PYF"

You can read more about using conda environments in the Managing Environments section of the conda documentation.  

### Steps to reproduce the main analysis in the publication
Typical results can be reproduced by executing the Jupyter Python notebooks  
The simulated examples are divided into folders according to the product name

1. The code for strain reconstruction is placed in the folder of strain reconstruction for each example  
**code file:**strain reconstruction.ipynb

2. The code for biosynthesis-growth relationship simulation is placed in the folder of biosyntehsis-growth relationship for each example  
**code file:**biosynthesis strain -FBA.ipynb, using FBA constraint  
**code file:**biosynthesis strain -kinetics.ipynb, using FBA and kinetic constraints  
**code file:**biosynthesis strain -thermodynamics.ipynb, using FBA and thermodynamic constraints  
**code file:**biosynthesis strain -kinetics and thermodynamics.ipynb, using FBA, kinetic and thermodynamic constraints

3. The code for MDF thermodynamic analysis for biosynthesis pathways is placed in the folder of MDF thermodynamic analysis for each example  
**code file:**MDF thermodynamic analysis for product biosynthesis.ipynb

4. The code for production prediction is placed in the folder of production prediction for each example  
**code file:**production prediction for product biosynthesis -FBA.ipynb, using FBA constraint  
**code file:**production prediction for product biosynthesis -kinetics.ipynb, using FBA and kinetic constraints  
**code file:**production prediction for product biosynthesis -thermodynamics.ipynb, using FBA and thermodynamic constraints  
**code file:**production prediction for product biosynthesis -kinetics and thermodynamics.ipynb, using FBA, kinetic and thermodynamic constraints

5. The code for metabolic engineering strategy optimization is placed in the folder of metabolic engineering strategy optimization for each example  
**code file:**metabolic engineering strategy optimization for product biosynthesis.ipynb  
