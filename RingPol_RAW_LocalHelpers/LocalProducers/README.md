# LocalProducers

**Local massive-multithreading TableProducer framework for the RingPol analysis.**

This folder lives inside:
```
PhD_codes/RingPol_RAW_LocalHelpers/LocalProducers/
```

It provides a **fully local alternative to the Hyperloop production workflow**, allowing you to:

- Test local production before escalating to Hyperloop (we need to make sure everything works before requesting pulls in GitHub!)
- Run ALICE O2 TableProducers locally with **extreme parallelism** (200–300+ threads)
- Efficiently stage data from a slow LAN to a fast local SSD
- Produce AO2Ds and QA outputs in a **Hyperloop-compatible format**
- Seamlessly plug the results into the legacy `DerivedDataHY` consumer framework

---

## Core idea

> Replace the distributed Hyperloop production stage with a **controlled, high-performance local pipeline**, while keeping the downstream analysis completely unchanged.

---

## 📋 Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Directory structure](#2-directory-structure)
3. [Scripts at a glance](#3-scripts-at-a-glance)
4. [Workflow A — Full local production](#4-workflow-a--full-local-production)
5. [Workflow B — Recovery and iteration](#5-workflow-b--recovery-and-iteration)
6. [Script reference](#6-script-reference)

---

## 1. Prerequisites

All scripts depend on the **O2/O2Physics environment**.

```bash
alienv enter O2Physics/latest
```

or pinned:

```bash
alienv enter O2Physics/v<TAG>-1
```

Verify:

```bash
which o2-analysis-lf-lambdajetpolarizationions
which hadd
which root
```

> **Important:**  
> If you open a new terminal mid-workflow, you must re-enter `alienv`.

---

## 2. Directory structure

Expected layout:

```text
<STORAGE_BASE>/
├── <DATASET_NAME>/
│   ├── input_data_storage.txt
│   ├── <WAGON_SHORTNAME>/
│   │   ├── AO2Ds/
│   │   └── results_producer/
│   │       ├── logs/
│   │       └── used_configs/
├── producer_configs/
└── consumer_configs/
```

> **Rule of thumb:**  
> Each wagon folder is **self-contained and reproducible**.

---

## 3. Scripts at a glance

| Script | You call it? | Purpose |
|-------|-------------|--------|
| `init_localWagon_structure.sh` | ✅ Yes | Auto-creates folder structure from configs |
| `run_all_producers.sh` | ✅ Yes | Runs all producer configurations |
| `runProducerAsHY.sh` | 🔧 Called internally | Core monolithic orchestrator |
| `merge_local_producer_qa.sh` | ⚠️ Fallback | Manual QA merging if needed |
| `run_local_consumers.sh` | ✅ Yes | Runs legacy consumer pipeline locally |

---

# 4. Workflow A — Full local production

Use this when starting from scratch.

---

### Step 1 — Scaffold the wagon structure

```bash
./init_localWagon_structure.sh
```

- Scans `producer_configs/`
- Infers wagon names automatically
- Creates all required directories
- Safe to re-run (uses `.initialized` sentinel files)

---

### Step 2 — Prepare input data

Ensure your dataset contains:

```text
input_data_storage.txt
```

This must list all input files (`file:...` paths).

---

### Step 3 — Run all producers

```bash
./run_all_producers.sh
```

This is the **main entry point**.

It will:

- Discover all `dpl-config-*.json`
- Loop over configurations
- Call `runProducerAsHY.sh` for each one
- Perform pre-flight checks
- Track progress, timing, and failures
- Handle interrupts safely (`Ctrl+C`)
- Print a final summary table
- Supports running on a fraction of the dataset via a percentage argument. Extremely useful for quick validation runs!

> ⚠️ **This step is extremely hardware-intensive**  
> Expect full CPU, RAM, and disk usage.

---

### Step 4 — Inspect outputs

Each wagon produces:

```text
AO2Ds/
results_producer/
```

If everything succeeded, you are ready for consumers.

---

### Step 5 — Run consumers

```bash
./run_local_consumers.sh
```

This will:

- Auto-discover wagon folders
- Generate a fake registry
- Call the legacy `DerivedDataHY` framework
- Run:
  - `extractDeltaErrors.cxx`
  - `signalExtractionRing.cxx`

---

# 5. Workflow B — Recovery and iteration

---

### B1 — QA merge failed

```bash
./merge_local_producer_qa.sh <WORK_DIR>
```

Use this if the final `hadd` step failed.

> This avoids re-running the entire producer pipeline.

---

### B2 — Re-run consumers only

```bash
./run_local_consumers.sh
```

No need to touch producers.

---

### B3 — Add new producer configs

1. Add new `dpl-config-*.json` to:
   ```text
   producer_configs/
   ```

2. Re-run:
   ```bash
   ./init_localWagon_structure.sh
   ./run_all_producers.sh
   ```

---

# 6. Script reference

---

## `init_localWagon_structure.sh`

```bash
./init_localWagon_structure.sh [OUTPUT_DATASET_DIR] [PRODUCER_CONFIGS_DIR]
```

**Purpose:**  
Automatically scaffolds wagon directories based on available configs.

**Key behavior:**

- Extracts wagon name from config filename
- Creates:
  ```text
  <WAGON>/
  ├── AO2Ds/
  └── results_producer/
      ├── logs/
      └── used_configs/
  ```
- Uses `.initialized` sentinel → safe to re-run

---

## `run_all_producers.sh`

```bash
./run_all_producers.sh [PERCENTAGE] [OUTPUT_DATASET_DIR] [INPUT_LIST] [PRODUCER_CONFIGS_DIR]
```

**Purpose:**  
The global wrapper for the production phase. It iterates over all producer configurations found in `PRODUCER_CONFIGS_DIR` and calls the `runProducerAsHY.sh` orchestrator for each one.

---

### **Arguments**

- **PERCENTAGE (Optional, default = 100)**  
  Fraction of the input dataset to process.  
  Accepts integers or floating-point values (e.g., `1`, `10`, `0.1`).  
  This is extremely useful for quick pipeline validation without processing the full dataset.

- **OUTPUT_DATASET_DIR (Optional)**  
  Target directory where the local Hyperloop-style structure will be created.

- **INPUT_LIST (Optional)**  
  Text file containing `file:` paths to process.

- **PRODUCER_CONFIGS_DIR (Optional)**  
  Directory containing your `dpl-config-*.json` files.

---

### **Dataset Sampling**

When `PERCENTAGE < 100`, the script:

1. Counts the number of valid entries in `INPUT_LIST`
2. Computes the corresponding fraction
3. Creates a temporary truncated input list
4. Uses this reduced list for all producers

The original input list is **never modified**.

At least **one file is always processed**, even for very small percentages.

---

### **Examples**

Run full dataset (default):
```bash
./run_all_producers.sh
```

Run only 1% of the data:
```bash
./run_all_producers.sh 1
```

Run 0.1% (very fast sanity test):
```bash
./run_all_producers.sh 0.1
```

Custom paths + sampling:
```bash
./run_all_producers.sh 5 /data /input.txt /configs
```

---

### **Additional Features**

- **Progress tracking**
  - Per-job runtime
  - Total elapsed time
  - Percentage of configs processed
  - Estimated time remaining (ETA)

- **Graceful interrupt handling**
  - `Ctrl+C` safely stops execution
  - Partial failure summary is printed

- **Failure tracking**
  - All failed configurations are collected
  - Final output includes a formatted failure table

- **Robust pre-flight checks**
  - Verifies required files and directories before execution

---

### **Notes**

- Sampling is **deterministic** (uses the first N entries of the list).  
  This ensures reproducibility across runs.

- If random sampling is desired, the script can be trivially modified to use `shuf`.

- The sampling mechanism is applied **once globally**, meaning all producer configurations operate on the exact same subset of data.

---

## `runProducerAsHY.sh`

```bash
./runProducerAsHY.sh <OUTPUT_DATASET_DIR> <INPUT_LIST> <CONFIG> <AOD_WRITER>
```

**Purpose:**  
The **core monolithic orchestrator**.

**Pipeline:**

1. Split input into batches (~50 files)
2. Stage data to SSD (LAN → local)
3. Run O2 TableProducer (200+ threads, large SHM)
4. Store resulting `AO2D_*.root`
5. Delete staged data
6. Merge QA outputs
7. Generate `input_data_storage.txt`

> Designed to saturate high-end hardware while avoiding I/O bottlenecks.

---

## `merge_local_producer_qa.sh`

```bash
./merge_local_producer_qa.sh <WORK_DIR>
```

**Purpose:**  
Manual fallback for QA merging.

**Behavior:**

- Merges:
  ```
  AnalysisResults_*.root → AnalysisResults_merged.root
  ```
- Preserves batch files if failure occurs

---

## `run_local_consumers.sh`

```bash
./run_local_consumers.sh [DATASET_DIR] [CONSUMER_CONFIGS_DIR]
```

**Purpose:**  
Bridge to legacy `DerivedDataHY`.

**Key mechanism:**

- Auto-discovers wagons
- Generates fake registry (`local_registry_auto.conf`)
- Calls:
  ```bash
  runDerivedDataConsumer_HY.sh
  ```
- Runs ROOT macros:
  - `extractDeltaErrors.cxx`
  - `signalExtractionRing.cxx`

> The legacy framework is completely unaware the data was produced locally.

---

## 🔑 Final insight

> This framework **decouples production from the Grid**, turning a distributed workflow into a **local, high-performance pipeline**, while preserving full compatibility with existing analysis tools. The idea is to have quick QAs before going to production level.

---