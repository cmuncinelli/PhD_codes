# Vorticity Hydro

## **Work in Progress**
> This folder contains proof-of-concept and development codes. Some files in `0Archive/` represent
> earlier approaches that were superseded, but are kept for reference.

Hydrodynamic and PYTHIA-level toy model codes for **Λ polarization from fluid vorticity**, used
to develop and validate the ring polarization observable before applying it to data.

The toy model provides a controlled environment where a known polarization signal can be injected,
allowing the ring observable machinery to be tested end-to-end.

## Main Analysis Codes

| File | Description |
|---|---|
| `lambda_pol_toy_model.cxx` | Main hydrodynamic toy model. Simulates Λ polarization from a vorticity field, reconstructs the ring observable, and validates the analysis chain. Supports parallelization via OpenMP. |
| `lambda_pol_toy_model_0_5cent.cxx` | Variant of the toy model restricted to the 0–5% centrality class. |
| `lambda_pol_PYTHIA.cxx` | Applies the ring polarization analysis to PYTHIA generator-level events. Since PYTHIA contains no polarization, this provides a zero-signal cross-check at the generator level. |
| `DebugGlobalToCompactIdx.cxx` | Debug utility for validating the global-to-compact index mapping used in the parallelization scheme. |
| `angular_distribution_sampler_test.py` | Python script to test and visualize the proton decay angular distribution sampler used in the toy model. |

## Reference and Notes

| File | Description |
|---|---|
| `run_lambda_pol_toy_model.txt` | Reference commands for running the toy model. |
| `run_lambda_pol_PYTHIA.txt` | Reference commands for running the PYTHIA polarization code. |
| `usual_input_angular_distribution_sampler.txt` | Typical inputs for the angular distribution sampler. |
| `Tracking_parallelization_performance.txt` | Notes on parallelization performance across development iterations. |
| `get_pythia_masses.txt` | Reference particle masses used in the PYTHIA-level code. |

## Archive

The [`0Archive/`](0Archive) folder contains earlier versions of the toy model that were superseded
but are kept for reference:

| File | Description |
|---|---|
| `lambda_pol_toy_model_TThreaded_Parallelization.cxx` | Parallelization attempt using ROOT's `TThreadedObject`. Superseded by the `std::vector` buffer approach. |
| `lambda_pol_toy_model-OpenMP_attempt.cxx` | Intermediate OpenMP parallelization version. |
| `lambda_pol_toy_model-ChatGPT_hallucination.cxx` | A test based on AI-generated code that failed to capture the analysis logic correctly. Kept as a reminder of the pitfalls of uncritical AI assistance. |