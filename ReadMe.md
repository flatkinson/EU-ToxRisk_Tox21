## Tox21 data visualization for EU-ToxRisk

This repo contains a series of Jupiter Notebooks to download the [Tox21](https://www.epa.gov/chemical-research/toxicology-testing-21st-century-tox21) data from PubChem (_via_ the [PUG API](https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html) and reformat it for visualisation in [DataWarrior](https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html).
       
DataWarrior is an open source tool similar to Spotfire and Vortex; for those new to DataWarrior, a comprehensive manual is available [here](http://www.openmolecules.org/help/basics.html). 

The notebooks are intended to be run in sequence; intermediate data structures are stored as pickle files in the directory `data/` and the final output file, `tox21_data_for_eu-toxrisk.csv`, is also written to this directory. This CSV file may be loaded into DataWarrior and the template `tox21_data_for_eu-toxrisk.dwat` (found the directory `files/`) applied to set up the desired visualization.

Briefly...

* `0_Tox21_assays.ipynb`: determines set of Tox21 Summary assays in PubChem
* `1_Tox21_compounds.ipynb`: retrieves set of substances and compounds tested in Tox21 assays from PubChem
* `2_Case_study_compounds.ipynb`: pre-processes list of EU-ToxRisk case study compounds
* `3_Tox21_compounds_vs_case_study_compounds.ipynb`: determines overlap of Tox21 compounds in PubChem and the case study compounds
* `4_Tox21_get_data.ipynb`: retrieves Summary assay data from PubChem
* `5_Tox21_reformat_data.ipynb`: reformats assay data for visualization in DataWarrior

Information on the EU-ToxRisk case strudy compounds is stored in the file `case_study_compounds_141216.txt` (found the directory `files/`); as the project progresses, this file may require updating. If this is the case, or if new Tox21 data is released, the visualization will be regenerated and the fact noted here.

The notebook `3_Tox21_compounds_vs_case_study_compounds.ipynb` requires a connection to a local instance of the ChEMBL database. If this is not available, the portion of the code using the connection should be commented out. If this is a problem, please let me know and I will modify the notebook to remove the need for database access by default.

Some of the notebooks can take a while to run, so, for convenience, the complete DataWarrior visualization file `tox21_data_for_eu-toxrisk.dwar` is available _via_ [DropBox](https://www.dropbox.com/sh/0iwk0ho59exhgkw/AABWLxsruT2huFpxvY-Ju77Pa?dl=0).

Note that this is a work in progress, so things may change. For example, it has been suggested taht some Tox21 assays may be of lesser utility than others and could thus be excluded or hidden by default. I have not yet looked into this, but will do so as time allows.

This project is based on preliminary work by Patricia Bento (patricia@ebi.ac.uk); the case study compound information was compiled by Anne Hersey (ahersey@ebi.ac.uk). Please address any comments or suggestions to me (francis@ebi.ac.uk).
