# ppigepred

Ppigepred (Protein-Protein Interaction GEne PREDiction) is a *very* simple software that takes a protein-pretein interaction network and a list of genes related to a phenotype of interest and predicts a list of other candidate genes that are potentially related to the phenotype.

Prediction is done using a basic random walk with restart approach.

# Installation

```bash
$ git clone https://github.com/mzdravkov/ppigepred
$ cd ppigepred

$ python3 -m venv env
$ source env/bin/activate

$ pip install -r requirements.txt
```

# Usage