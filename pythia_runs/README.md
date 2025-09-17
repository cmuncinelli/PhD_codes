A simple folder for running PYTHIA8 and PYTHIA Angantyr.
Mostly complimentary to the vorticity_hydro analysis in a sense that PYTHIA's ring observable is expected to be zero (all decays are isotropic and there is no explicit consideration of polarization in it)

- [**input_cards_ALICE**](input_cards_ALICE): All ALICE input cards come from the O2DPG repository (O2DPG/MC/config/ALICE3/pythia8/generator)

The [**PythiaGenMin.cxx**](PythiaGenMin.cxx): code is a direct adaptation of the PythiaGen.cxx code in the laguerre_project repository. It is just a parallel adaptation of that code, stripped of most of the unnecessary stuff that we required for comparison with the [1807.11321 [nucl-ex], pp 7 TeV](https://arxiv.org/abs/1807.11321) paper.