# RingPol RAW Local Helpers

Local analysis framework for the **ring polarization observable**, operating on both raw AOD data
and Hyperloop-derived data outside of the Hyperloop environment. This folder collects the tools
needed to download, process, and analyze data locally -- complementing the O2Physics tasks that run
inside Hyperloop.

## Subfolders

### [`Local_framework/`](Local_framework)
Scripts and configuration for running the O2 DPL analysis workflow locally on raw AOD files.
Includes download helpers, path generators, and multi-threaded execution scripts.

### [`DerivedDataHY/`](DerivedDataHY)
Framework for downloading and analyzing Hyperloop-produced derived data locally.
See [`DerivedDataHY/README.md`](DerivedDataHY/README.md) for full documentation.

### [`JsonExamples/`](JsonExamples)
Example DPL workflow configuration files (`.json`) for various analysis scenarios: both-hyperon
mode, Lambda-only, AntiLambda-only, and full QA with permissive pT cuts.

## C++ Analysis Codes

| File | Description |
|---|---|
| `signalExtractionRing.cxx` | Main signal extraction code for the ring polarization observable. Extracts the ring signal from invariant mass spectra across the (pT, centrality, Δφ) analysis bins. |
| `signalExtractionRingTest.cxx` | Test and development version of the signal extraction code. |
| `extractDeltaErrors.cxx` | Extracts and propagates statistical uncertainties using the Delta Method, which more accurately preserves the covariance terms to avoid SEM mathematical pitfalls. (Did not improve by much though) |
| `makeCumulativeDCAdauProfile.cxx` | Evaluates the robustness of the ring observable against varying Distance of Closest Approach (DCA) cuts, specially the DCA of the daughter particles wrt the PV. Condenses 2D/3D differential QA profiles into unified cumulative distributions. |
| `auxiliaryPlots.cxx` | Cross-configuration aggregation macro. Combines baseline, MC reference, Helicity Toy Model inputs, and systematic variations (e.g., Data-like Jet, Rand Jet) of ring observables into appropriately scaled, unified canvases for easy comparison. |
| `testVersionOfHeaders.h` | Development version of shared header definitions. |

## Utilities

| File | Description |
|---|---|
| `mergeAODDerived.sh` | Merges AOD derived data files prior to local analysis. |