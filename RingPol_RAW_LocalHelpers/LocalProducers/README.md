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
7. [Hardware tuning reference](#7-hardware-tuning-reference)

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
- Call `runProducerAsHY.sh` for each one (pinned to NUMA node 0 via `numactl`)
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

### **NUMA Pinning**

Each invocation of `runProducerAsHY.sh` is launched with:

```bash
numactl --cpunodebind=0 --preferred=0
```

This pins all processing to **NUMA node 0**, preventing cross-socket memory traffic. This is intentional and coordinated with the batch size tuning inside `runProducerAsHY.sh` (see [Hardware tuning reference](#7-hardware-tuning-reference)).

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
  - Estimated time remaining (ETA) computed as `avg_time_per_job × remaining_jobs`

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
./runProducerAsHY.sh [--TOF | --no-TOF] [--do-log] <OUTPUT_DATASET_DIR> <INPUT_LIST> <CONFIG_JSON> <AOD_WRITER_JSON>
```

**Purpose:**  
The **core monolithic orchestrator**.

**Pipeline:**

1. Split input into batches (20 files each — see [Hardware tuning reference](#7-hardware-tuning-reference))
2. Stage data to SSD (LAN → local), with an inline progress bar
3. Run O2 TableProducer with tuned parallelism and large SHM
4. Store resulting `AO2D_*.root`
5. Delete staged data (batch workspace freed immediately after each batch)
6. Merge QA outputs via `hadd -j 16`
7. Generate `input_data_storage.txt`

---

### **Flags**

#### `--TOF` / `--no-TOF` — TOF pipeline selection

Controls whether the TOF-extended analysis pipeline is used.

| Flag | Behavior |
|------|----------|
| *(omitted)* | **Auto-detect** (default): reads `"processDataWithTOF"` from the JSON config |
| `--TOF` | Force TOF pipeline on, regardless of config |
| `--no-TOF` | Force TOF pipeline off, regardless of config |

**When TOF is enabled**, the pipeline additionally includes:
- `o2-analysis-ft0-corrected-table`
- `o2-analysis-pid-tof-base` (with `--pipeline tof-event-time`)
- `o2-analysis-lf-strangenesstofpid` (with `--pipeline strangenesstofpid`)

**When TOF is disabled**, those three tasks are omitted entirely from the DPL pipeline.

The source of the TOF decision (auto-detected from config vs. forced via CLI) is printed at startup:
```
TOF Pipeline: TRUE (Source: JSON Config)
TOF Pipeline: FALSE (Source: Manual CLI Flag)
```

#### `--do-log` — Per-batch log files

| Flag | Behavior |
|------|----------|
| *(omitted)* | **Default**: all O2 output is redirected to `/dev/null` |
| `--do-log` | Saves full stdout+stderr of each batch to `results_producer/logs/batch_<ID>.log` |

Log files are useful for debugging pipeline failures but produce large outputs on full runs. The default is intentionally silent to reduce I/O overhead during production.

> **Debugging tip:** Add `--fairmq-rate-logging 1` to the last O2 task in the pipeline (in `runProducerAsHY.sh`) to diagnose DPL message-rate bottlenecks. See the comment in the script.

---

### **Resource configuration**

At startup the script auto-detects hardware and prints a resource summary, for example:

```
Resource configuration:
  Total cores        : 168
  Total threads      : 336
  Max usable cores   : 84
  Max usable threads : 168
  Used threads       : 42
  Readers            : 8
  I/O Threads        : 8
```

All tuning knobs live in the `TUNING KNOBS` block near the top of the script and are documented in [Section 7](#7-hardware-tuning-reference).

---

### **Staging progress bar**

During the stage-in phase (LAN → SSD), an inline ASCII progress bar is displayed per batch:

```
Progress: [######################        ]  73% (15/20)
```

---

### **Batch workspace isolation**

Each batch runs in its own subdirectory under `WAGON_DIR/temp_staging_area/work_<ID>/`. The workspace is **deleted immediately** after the batch completes, keeping the SSD footprint at one batch at a time.

The entire `temp_staging_area/` is also cleaned on exit (normal or interrupted) via a `trap ... EXIT` handler.

---

### **Final merge**

QA files are merged with:

```bash
hadd -f -j 16 AnalysisResults_merged.root @merge_list.txt
```

The `-j 16` flag enables 16 parallel merge threads. Individual per-batch `AnalysisResults_*.root` files are deleted after a successful merge.

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

# 7. Hardware tuning reference

This section documents all tunable knobs in `runProducerAsHY.sh` and the rationale behind the current values. All are defined in the `TUNING KNOBS` block near the top of the script.

---

### Hardware target

The current defaults are tuned for:

- **CPU:** 2× AMD EPYC 9634 (84 cores / 168 threads per socket)
- **RAM:** 512 GB total
- **SSD:** Micron 7450 MTFDKBG3T8TFR, ~5000 MB/s read

---

### Batch size: `FILES_PER_BATCH=20`

Each batch consists of **20 input files**, corresponding to roughly 60 GB of raw data.

**Rationale (NUMA-aware):**  
Since `run_all_producers.sh` pins execution to NUMA node 0 (`numactl --cpunodebind=0 --preferred=0`), only ~257 GB of the 512 GB RAM is available to the process. 20 files at ~3 GB each keeps the entire LAN-to-SSD transfer resident in the Linux Page Cache within node 0's memory budget, without competing with the 128 GB SHM allocation or forcing cross-socket memory spills.

> For quick debugging runs, 10 files per batch is also acceptable.

---

### Half-core policy

```bash
MAX_CORES_ALLOWED=$((TOTAL_CORES / 2))
MAX_THREADS_ALLOWED=$((MAX_CORES_ALLOWED * 2))
```

Only half of the physical cores are used. This is intentional for a **shared server** context: it prevents monopolizing the machine while still saturating memory bandwidth.

---

### Pipeline parallelism

The DPL pipeline is split into independent processes, each with a configurable worker count (`--pipeline <task>:<N>`). Current values:

| Task | Workers | Notes |
|------|---------|-------|
| `eventselection-run3` | 4 | Lightweight flag-checking |
| `mult-cent-table` | 4 | Was a bottleneck; keep moderate |
| `propagation-service` | 8 | Propagation workload |
| `pid-tpc-service` | 10 | CPU bottleneck; highest priority |
| `tof-event-time` | 4 | TOF pipeline only |
| `strangenesstofpid` | 4 | TOF pipeline only |
| `lambdajetpolarizationions` | 8 | Main analysis task |
| **Total** | **42** | (without TOF: 34) |

**Why reduced from prior iteration?**  
Profiling showed `sys > user` time, indicating the OS was overwhelmed by ZMQ inter-process communication overhead from too many concurrent DPL workers. The process count was reduced significantly to bring IPC overhead under control.

> The commented-out block in the script (`PIPE_EVENTSEL=10`, `PIPE_PROPAGATION=22`, ...) preserves the previous values for reference.

---

### I/O tuning

| Parameter | Value | Notes |
|-----------|-------|-------|
| `READERS` | 8 | Number of AOD reader processes |
| `IO_THREADS` | 8 | Writer-side I/O threads |
| `SHM_SIZE` | 128 GB | Shared memory segment for DPL messaging |
| `MEM_RATE_LIMIT` | 8 GB/s | AOD memory rate limit, scaled for 20-file batches |

**Reader safeguard:**  
If `REQUESTED_READERS > FILES_PER_BATCH`, readers are automatically capped to `FILES_PER_BATCH`. This prevents a crash caused by requesting more reader processes than there are input files in a batch. A warning is printed:

```
[I/O Tuning] WARNING: Requested readers (8) exceeds files per batch (5). Capping READERS to 5.
```

> **Spawners** (`--spawners`) were considered but disabled. They appear to fail when configuring DPL links in heavily pipelined workflows. The relevant commented-out line is kept in the script for future reference.

---

## 🔑 Final insight

> This framework **decouples production from the Grid**, turning a distributed workflow into a **local, high-performance pipeline**, while preserving full compatibility with existing analysis tools. The idea is to have quick QAs before going to production level (i.e., before doing PRs, but still having a good amount of statistics locally).

---