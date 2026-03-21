# PYTHIA Runs

**WARNING: Work in Progress.**
Code and results in this folder are under active development. Some optimization studies and
profiling results may reflect intermediate stages of the code.

--

A folder for running PYTHIA8 and PYTHIA Angantyr simulations. Mostly complementary to the
[vorticity_hydro](../vorticity_hydro) analysis in the sense that PYTHIA's ring observable is
expected to be zero -- all decays are isotropic and there is no explicit consideration of
polarization. Since PYTHIA does not include spin-orbit coupling or polarization dynamics by
construction, any non-zero signal measured here would indicate a bias or systematic effect in the
analysis framework.

FastJet is used for jet finding within the generator-level events.

--

## Main Generator Code

- [**PythiaGenMin.cxx**](PythiaGenMin.cxx): Direct adaptation of the `PythiaGen.cxx` code from
  the laguerre_project repository (personal repository for early undergraduate-level studies I performed). A parallelized (OpenMP) version stripped of most of the
  overhead not needed for this analysis, producing minimal-footprint output trees with only the
  particle information relevant to the ring polarization analysis (Lambda, LambdaBar, p, pBar,
  jet kinematics). (Removed the unnecessary stuff that we required for comparison with the [1807.11321 [nucl-ex], pp 7 TeV](https://arxiv.org/abs/1807.11321) paper on the undergrad research code)
- [**CentralityCorrection_FixingOverflows.cxx**](CentralityCorrection_FixingOverflows.cxx):
  Post-processing correction for centrality estimation overflow bins.
- [**Makefile**](Makefile): Build system for the generator codes.

## Input Cards

- [**input_cards_ALICE/**](input_cards_ALICE): PYTHIA8 configuration cards for the collision
  systems of interest. All ALICE input cards originate from the O2DPG repository
  (O2DPG/MC/config/ALICE3/pythia8/generator).

| File | System |
|---|---|
| `pythia8_OO.cfg` | O--O collisions |
| `pythia8_PbPb_536tev.cfg` | Pb--Pb at sqrt(sNN) = 5.36 TeV |
| `pythia8_PbPb_536tev_AllowLambdaDecays.cfg` | Same, with Lambda decays enabled |
| `pythia8_pp_136tev.cfg` | pp at sqrt(s) = 13.6 TeV |
| `pythia8_pp_136tev_AllowLambdaDecays.cfg` | Same, with Lambda decays enabled |

## Profiling and Optimization

The generator was iteratively optimized over several passes (see commit history). Performance
studies are collected in [profiling/](profiling) and [runtime_studies/](runtime_studies).

| File | Description |
|---|---|
| `profiling/profile_perf.sh` | Profiles the binary with Linux `perf`. |
| `profiling/extract_top_functions_from_profile_perf.sh` | Summarizes the hottest functions from a `perf` report. |
| `profiling/run-callgrind.sh` | Runs Valgrind Callgrind for detailed call graph profiling. |
| `profiling/simpler_perf.sh` | Lightweight perf wrapper for quick profiling passes. |
| `pgo_build_and_run.sh` | Profile-guided optimization (PGO) build script (experimental). |
| `Choosing_compiler_optimizations.txt` | Notes on compiler flag choices and their measured effects. |
| `testing_optimizations.txt` | Raw notes from optimization testing sessions. |
| `runtime_studies/` | Collected runtime measurements across optimization passes. |

## Experimental / Archive

| File | Description |
|---|---|
| `PythiaGenMin-WeirdTestWithPythiaParallel.cxx` | Experimental attempt at using `PythiaParallel` (superseded by the OpenMP approach). |
| `PythiaGenMin-StructTestForMemoryContinuity.cxx` | Memory layout experiment using structs for cache-friendlier access. |

## Typical Run Commands

See `usual_input_PythiaGen.txt` and `usual_input_CentralityCorrection_FixingOverflows.txt` for
reference invocation examples.