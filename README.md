# PhD Codes -- Cicero Muncinelli

## **Work in Progress -- Please Read**
> This repository is an active, evolving research codebase spanning the full duration of a PhD (2025--2030).
> Although every effort is made to keep the code correct and well-documented, some parts may be
> unfinished, outdated, or in the process of being reworked. Earlier versions of some codes have been
> deliberately preserved for bookkeeping purposes -- they may contain naive approaches or undeveloped
> ideas that could still prove useful later. No guarantee of correctness (or more precisely, of *completeness*) is implied
> for any file in this repository.

---

## Overview

This repository collects most of the analysis and simulation codes developed during my PhD at
[IFGW -- Instituto de Física Gleb Wataghin](https://portal.ifi.unicamp.br/), Universidade Estadual
de Campinas (Unicamp), within the ALICE collaboration at CERN.

The central physics goal is to **search for local Λ (and anti-Λ) polarization induced by jets** -- in a search for Hydro-predicted vorticity ring structures that may form around them -- **in heavy-ion and hadronic
collisions** (Pb--Pb, O--O, and pp). PYTHIA and PYTHIA Angantyr Monte Carlo simulations, which
contain no polarization signal by construction, serve as essential cross-checks. Some codes in this
repository are more general-purpose and may not be directly tied to this specific analysis.

The O2Physics tasks that run inside the ALICE Hyperloop framework live in a separate repository and
are not included here.

(For those, look at my fork at [cmuncinelli/O2Physics](https://github.com/cmuncinelli/O2Physics)!)

---

## Repository Structure

| Folder | Description |
|---|---|
| [`vorticity_hydro/`](vorticity_hydro) | Hydrodynamic and PYTHIA-based toy models for Λ polarization from fluid vorticity. Used as proof-of-concept and ring observable development. |
| [`pythia_runs/`](pythia_runs) | PYTHIA8 and PYTHIA Angantyr generator codes with FastJet integration. Produces zero-polarization MC reference samples for the ring observable analysis. |
| [`RingPol_RAW_LocalHelpers/`](RingPol_RAW_LocalHelpers) | Local analysis framework for the ring polarization observable. Includes tools for downloading, processing, and analyzing both raw AOD and Hyperloop derived data. |
| [`LambdaV0Radius_QA/`](LambdaV0Radius_QA) | Quality assurance study for systematic effects related to the V0 decay radius in Λ reconstruction. |
| [`IC_helper_tasks/`](IC_helper_tasks) | Utility scripts written to support undergraduate research students. |

---

## Notes on Paths and Environment

Some scripts contain hardcoded paths such as `/home/cicero/` or references
to cluster nodes (e.g. `jarvis*`). These reflect the original computing
environment at Unicamp and should be adapted to your local setup.

---

## Acknowledgments

This work is supported by the Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP).

I would like to acknowledge the people of the state of São Paulo, whose public investment in science makes research like this possible -- even when its impact is not immediately visible. I am sincerely grateful for that support and trust.

This research is carried out at the Instituto de Física Gleb Wataghin (IFGW), at the Universidade Estadual de Campinas (Unicamp).