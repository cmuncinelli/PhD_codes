A simple collection of tools that may come in hand before creating a PR (see [PRs/](PRs/)) or for shared memory cleanup if O2Physics is being stubborn with memory management (that shouldn't be the case in newer versions of the framework, nor of my own bash scripts that call it and now wait for graceful exits after Ctrl+C is issued).

Contents:

- [jsonCleaner.py](jsonCleaner.py) -- prunes a DPL configuration JSON down to the values that are genuinely non-default.
- [shrMemCleanup.md](shrMemCleanup.md) -- commands for clearing stuck shared memory after an ungraceful exit.
- [PRs/](PRs/) -- linting, clang-format and cppcheck helpers to run before opening a PR.

---

# jsonCleaner.py

> Be **VERY CAREFUL** when running this tool! I cannot stress this enough, honestly. If you change your defaults afterwards, you can silently corrupt a bunch of JSON configurables. Say that, for instance, `Configurable<bool> analyseLambda` was set to true in your defaults. All JSON keys that said `true` for `JustLambda` workflows were deleted after you ran this. If you then change the default to `false` and re-compile O2Physics, your `JustLambda` workflows will all be thrown away!

> This file is not meant to be ran every single time you update your O2Physics code. It is meant as a clean starting point whenever you want to build your first clean JSON from the automatically generated scaffold. Be VERY CAREFUL when running it in already-established workflows!!!

The json cleaner is a new tool that allows you to prune the configurables laid out in a JSON after you create a template for it via David's trick for JSON automatic generation. See slide 7 of [David's DDChinellato-O2AT5-HandsOn-01.pdf hands-on session](https://indico.cern.ch/event/1574136/contributions/6785565/attachments/3170125/5636806/DDChinellato-O2AT5-HandsOn-01.pdf) for what I mean.

Essentially, it removes all default parameters from the `.json` file and leaves only the non-defaults that you've set by hand when running your consumer/producer locally. It reads the defaults straight out of the task's `.cxx`, so the reference is always the code that is actually compiled, never a transcription of it.

It also removes keys that have **no counterpart in the source at all**. Those matter more than they look: DPL only ever queries the options a device registered, so a stale key is inert, but a **renamed** one is worse than inert -- the old name is ignored and the new name quietly falls back to its compiled default.

## Requirements

- `python3` (standard library only, no pip installs). Tested on 3.12.
- **o2env is only needed when a default is not a plain literal** -- see [Symbol resolution](#symbol-resolution) below. Routine runs work outside it.

## Usage

Dry run (the default -- nothing is ever written without `--apply`):

```bash
python3 jsonCleaner.py \
  --source ~/alice/O2Physics/PWGLF/Tasks/Strangeness/lambdaJetPolarizationIonsDerived.cxx \
  --configs ~/RingPol/consumer_configs/
```

Same thing, actually applied to the JSON:

```bash
python3 jsonCleaner.py \
  --source ~/alice/O2Physics/PWGLF/Tasks/Strangeness/lambdaJetPolarizationIonsDerived.cxx \
  --configs ~/RingPol/consumer_configs/ \
  --apply
```

`--configs` takes a folder or a single `.json`. Backups of every file that actually changes go to `<configs>/prePruningJsons/<YYYYmmdd-HHMMSS>/`, one folder per run, so a second run never clobbers the first backup. The tool skips its own backup folder, of course.

> JSONs are small enough that keeping previous versions of them is not a significant cost in storage. In case this ever bothers you, just delete the older backups!

Writing the full compiled defaults out as a reference:

```bash
python3 jsonCleaner.py --source ~/alice/O2Physics/PWGLF/Tasks/Strangeness/lambdaJetPolarizationIonsDerived.cxx \
  --dump-defaults ~/PhD_codes/RingPol_RAW_LocalHelpers/JsonExamples/ --defaults-label DerivedConsumer
python3 jsonCleaner.py --source ~/alice/O2Physics/PWGLF/TableProducer/Strangeness/lambdaJetPolarizationIons.cxx \
  --dump-defaults ~/PhD_codes/RingPol_RAW_LocalHelpers/JsonExamples/ --defaults-label TableProducer
```

That produces `dpl-config-DerivedConsumer-Defaults.json` and `dpl-config-TableProducer-Defaults.json`. These are **reference dumps, not runnable configs**: they contain only the task's own device block, without the `internal-dpl-*` scaffolding a real workflow needs.

## What the report tells you

| Verdict | Meaning | Effect |
| --- | --- | --- |
| `KEEP` | genuinely differs from the compiled default | stays in the file |
| `PRUNE` | equal to the compiled default | removed |
| `NEAR` | differs, but by less than the safety net | removed, **with a warning** |
| `ORPHAN` | no such configurable in the source | removed |

Every device key that does not belong to the analysed source is listed under `untouched` and copied over byte for byte. This is what makes the tool safe on a producer config like `dpl-config-ITSandTPC.json`, where `lambdajetpolarizationions` sits alongside eleven other devices (`eventselection-run3`, `mult-cent-table`, `strangenesstofpid`, ...). A JSON containing no device from the given source is warned about and skipped, not modified.

## How the defaults are read

Four patterns are recognised inside any top-level `struct` in the source:

- `Configurable<T> name{"key", default, "help"};`
- `ConfigurableAxis name{"key", {...}, "help"};` (`VARIABLE_WIDTH` is the plain `0` sentinel of the stored vector)
- `struct : ConfigurableGroup { std::string prefix = "groupName"; ... }` -- the **`prefix` member**, not the instance name, is what becomes the JSON sublevel
- `PROCESS_SWITCH(struct, name, "help", default)` and its `_FULL` variant -- easy to forget, but a process switch is a configurable like any other

Declarations spanning several lines are handled, and comments are stripped in a literal-aware pass, which is not optional: `ccdbUrl`'s default is `"http://alice-ccdb.cern.ch"`, and a naive stripper truncates it at the `//`. Commented-out configurables are correctly ignored -- there are plenty in both of my sources, and they are one of the ways an orphan key gets created in the first place.

The device key is derived from the struct name with O2's task-name rule (lowercase, dash before each inner capital), and an explicit `TaskName{"..."}` in `adaptAnalysisTask` overrides it. Several structs in one file give several device keys, all handled in the same pass.

### Type handling

Defaults are canonicalized the way C++ would actually store them, with **two casts**, because narrowing happens in both directions in real code:

```cpp
Configurable<double> radiusJet{"radiusJet", 0.4f, ...};  // float literal widened to double
Configurable<float>  v0cospa  {"v0cospa",   0.995, ...};  // double literal narrowed to float
```

So the literal is first evaluated in its own type (`f` suffix means `float`), then cast to the declared type. `ConfigurableAxis` is a `Configurable<std::vector<double>>`, so only the first cast bites there -- which is why `0.1f` lands in the JSON as `0.10000000149011612`.

<a name="symbol-resolution"></a>
### Symbol resolution (the one thing that needs o2env)

Some defaults are not literals: `constants::math::PI`, `constants::math::TwoPI / 18`, `kCentFT0M`, `kAntiKt`. These are **not** hardcoded here, on purpose -- `constants::math::PI` is a `float`, so writing `math.pi` in Python would be wrong in the 8th digit and would flag every PI-based axis as a non-default.

Instead, the whole expression (not just the symbol) is handed to ROOT, in a generated probe macro that includes `CommonConstants/MathConstants.h`, `CommonConstants/PhysicsConstants.h` and the `enum` blocks lifted verbatim from the source being analysed. ROOT does the enum numbering and the arithmetic; the cast to `double` happens only at the print boundary, so the expression's own type is preserved.

This requires being inside o2env:

```bash
alienv enter O2Physics/latest
```

Results are cached in `symbolCache.json` next to the script, so **o2env is only needed the first time a new symbol appears**. If a symbol is missing and `root` is not on the PATH, the tool stops and says so rather than guessing. Point `--symbol-cache` elsewhere if you want a per-project cache.

## Numerical comparison, and why it has three bands

The DPL JSON writer serializes with `max_digits10` **of the declared type**, so it never round-trips exactly:

- a `double` loses about one ulp -- `1.2f` as an axis edge comes back as `1.200000047683716` instead of `1.2000000476837158` (relative ~2e-16)
- a `float` is cut to 9 significant digits -- `etaCut` is written as `0.899999976` (relative ~5e-9)

Bit-exact comparison is therefore impossible, and a single tolerance would either flag every float as non-default or hide real changes. The comparison uses three bands instead:

| Relative difference | Verdict |
| --- | --- |
| below the writer's noise floor (1e-8 for `float`, 1e-12 for `double` and axes) | `PRUNE`, silently |
| between that and **1e-6** | `NEAR`: **pruned anyway, with a warning** |
| above 1e-6 | `KEEP` |

Integers, booleans and strings are compared exactly. A default of exactly `0` also requires exact equality, since a relative tolerance means nothing against zero.

**The `NEAR` policy is deliberate: the key is removed, and it is on you to put it back.** The warning prints both values and the relative difference, on stderr, in dry run as well as under `--apply`. If the difference was intentional, restore that one key from `prePruningJsons/`. The reasoning is that a sub-ppm difference in a cut value or a bin edge is never physics, so the common case should be a clean file, with the rare case made loud rather than silent. `--eps-warn` moves the 1e-6 boundary if you disagree on a given run.

This band is what separates truncation noise from real discrepancies. On my producer config it silences nine harmless 9-digit truncations and leaves exactly one genuine warning: `jetConfigurations.radiusJet`, where the JSON holds a plain `0.4` while the source now says `0.4f` in a `Configurable<double>`.

## The safety invariant

Before anything is written -- in dry run too -- the tool checks, for every device and every configurable, that

```
merge(defaults, pruned)  ==  merge(defaults, original)
```

that is, that the effective configuration the device will see is unchanged, except for `NEAR` keys, which are allowed to move within the warning band. It also re-checks that every key called an orphan really is absent from the source. Any violation is fatal: it prints `FATAL:`, writes nothing, and exits with status 1, so a wrapper script can trap it.

This is not a flag and cannot be turned off. A tool whose failure mode is silently deleting a setting you meant to keep does not get to skip its own check.

## Known limits

- One source file per invocation. Run it once for the producer and once for the consumer.
- Only `Configurable`, `ConfigurableAxis`, `ConfigurableGroup` and `PROCESS_SWITCH` are recognised. An unsupported `Configurable<T>` is fatal rather than skipped -- better a stop than a wrong prune.
- Configurables inherited from an included header are not followed. Neither of my two tasks has any; if that changes, the tool will report the keys as orphans, which is the loud failure and not the quiet one.
- The `internal-dpl-aod-reader` block also carries framework defaults, but it belongs to no analysed source and is left alone by design.