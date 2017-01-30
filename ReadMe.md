## Tox21 data visualization for EU-ToxRisk

This repo contains a series of Jupiter Notebooks to download the [Tox21](https://www.epa.gov/chemical-research/toxicology-testing-21st-century-tox21) data from PubChem (_via_ the [PUG API](https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html) and reformat it for visualisation in [DataWarrior](https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html).
       
DataWarrior is an open source tool similar to Spotfire and Vortex; for those new to DataWarrior, a comprehensive manual is available [here](http://www.openmolecules.org/help/basics.html). 

The notebooks are intended to be run in sequence; intermediate data structures are stored as pickle files in the directory `data/` and the final output file, `tox21_data_for_eu-toxrisk.csv`, is also written to this directory. This CSV file may be loaded into DataWarrior and the template `tox21_data_for_eu-toxrisk.dwat` (found the directory `files/`) applied to set up the desired visualization.

Some of the notebooks can take a while to run, so, for convenience, the DataWarrior visualization is available _via_ [DropBox](https://www.dropbox.com/sh/0iwk0ho59exhgkw/AABWLxsruT2huFpxvY-Ju77Pa?dl=0).

Information on the EU-ToxRisk case strudy compounds is stored in the file `case_study_compounds_141216.txt` (found the directory `files/`); as the project progresses, this file may require updating. If this is the case, or if new Tox21 data is release, the visualization will be regenerated and the fact noted here.

Note that the notebooks are based on preliminary work by Patricia Bento (patricia@ebi.ac.uk); the case study compound information was compiled by Anne Hersey (ahersey@ebi.ac.uk). Please address any comments or suggestions to me (francis@ebi.ac.uk).
