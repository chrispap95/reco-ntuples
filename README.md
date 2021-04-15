# Coffea related files for the HGCAL reco-ntuples. 

Since the branch `coffea-cmssw`, the format of the n-tuples are now
[coffea][coffea] compliant. This branch stores the various information required
to assist with working with coffea and related tools.

For the documentation for setting up a coffea environment, see the official
[coffea introduction documentation][coffea-intro]. For complete beginners to
coffea, see the [unofficial introduction to coffea][coffea-umd]

## List of files

- [`hgcreco_scheme.py`](hgcreco_schema.py): the python file containing the schema
  for reading out the ntuples using awkward arrays.
- [`hgcreco_explore.ipynb`](hgcreco_explore.ipynb): Example notebook for
  exploring the ntuple in coffea/awkward notation.

[coffea]: https://coffeateam.github.io/coffea/
[coffea-intro]: https://coffeateam.github.io/coffea/installation.html
[coffea-umd]: https://umdcms.github.io/CoffeaTutorial/
