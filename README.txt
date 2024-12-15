PYF
The scripts for production simulation and  strategy optimazation are in a Pyomo modelling framework.

Environment
The scripts were written and tested with Python 3.9. 
The core libraries essential for the pipeline including: Cobrapy toolkit: version --0.25.0； Pyomo package: version --6.4.0； Gurobi solver: version --10.0.1.

Software
The packages used to run the pipeline was listed in requirements.txt. 
To install the requirements using pip, run the following code at command-line: $ pip install -r requirements.txt. 
To create a stand-alone environment named PYF with Python 3.9 and all the required package versions (especially for cobrapy is also available), run the following code: 
$ conda create -n PYF python=3.9 
$ conda activate PYF 
$ pip install -r requirements.txt 
$ python -m ipykernel install --user --name PYF --display-name "PYF" 
You can read more about using conda environments in the Managing Environments section of the conda documentation.

Steps to reproduce the main analysis in the publication
Typical results can be reproduced by executing the Jupyter Python notebooks:
The -FBA.ipynb, -kinetics.ipynb, -thermodynamics.ipynb, and -kinetics and thermodynamics.ipynb files in PYF folder simulate the biosynthesis pathway expression degrees.
The thermodynamic analysis for product biosynthesis.ipynb files in PYF folder simulate the relationships between biosynthesis tasks and growth interests. 
The two strain consortium for product biosynthesis -FBA.ipynb, two strain consortium for product biosynthesis -kinetics.ipynb, two strain consortium for product biosynthesis -thermodynamics.ipynb, and two strain consortium for product biosynthesis -kinetics and thermodynamics.ipynb files in PYF folder simulate the productions.
The strategy analysis for product biosynthesis.ipynb files in PYF folder optimize the biosynthesis strategies.
Download the ETGEMs_function.py at https://github.com/tibbdc/ETGEMs, and place the file in the codes folder. Rename the ETGEMs_function.py to ETGEMs_function_ETG.py