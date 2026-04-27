# Workflow Specification

This document defines the **standard stage-gated workflow** for high-throughput luminescent materials screening in this repository.

It is written for both:
- **agents** that need an explicit execution contract
- **human operators** who want a compact view of what each stage is allowed to do

---

## 1. Goal

The workflow is designed to:
1. generate donor/acceptor candidates under explicit topology rules
2. remove structurally unstable molecules early
3. use a **cheap proxy photophysics stage** to narrow the search space
4. reserve expensive TDDFT validation for only the best surviving candidates

This workflow supports:
- **TADF**
- **Fluorescence**
- **Phosphorescence**

---

## 2. Core design principle

The workflow is **stage-gated**.

That means:
- each stage has a clear purpose
- each stage has defined inputs and outputs
- a later stage must not silently compensate for a failed earlier stage
- if a gate fails, the workflow must stop and report the reason

In short:

`initialize -> generate -> stabilize -> proxy filter -> TDDFT validate -> summarize`

---

## 3. Stage overview

| Stage | Name | Purpose | Typical Cost | Main Output |
|------|------|---------|--------------|-------------|
| 0 | Initialization | define boundary conditions | very low | workflow config |
| 1 | Structure generation + stability | build candidates and remove unstable ones | low | stable xyz set |
| 2 | Semi-empirical photophysics proxy | estimate rough excitation/emission behavior | medium | shortlist |
| 3 | Ab initio validation | confirm elite candidates with TDDFT | high | final validated results |
| 4 | Spectrum + reporting | render outputs and summarize decisions | low | plots / tables / reports |

---

## 4. Stage 0 — Initialization (mandatory)

Run `scripts/screening_workflow_initializer.py` first.

### Responsibilities
- detect available compute environment
- assign `compute_tier`
  - `local_basic`
  - `local_gpu`
  - `remote_cluster`
- capture target constraints:
  - emission type
  - target wavelength window
  - spectrum-width preference
  - empirical Stokes shift
- discover available quantum engines
  - `gaussian`
  - `orca`
  - `qchem`
  - `pyscf`
- define runtime mode and missing-field inquiry rules

### Required output
A machine-readable workflow configuration that downstream stages can consume.

### Hard rule
Do **not** treat a ground-state vertical excitation as the final emission answer.

Final validation must follow:

`S0 optimization -> excited-state optimization (S1/T1 as appropriate) -> vertical emission`

---

## 5. Stage 1 — Structure generation and stability pre-screen

### Purpose
Generate chemically meaningful candidates and discard unstable structures early.

### Typical steps
1. assemble candidates with topology rules
2. generate 3D geometry with RDKit
3. run MMFF94 or equivalent cheap pre-optimization
4. export `.xyz`
5. run `xtb --opt` or equivalent geometry optimization

### Minimum acceptance criteria
A candidate must:
- assemble successfully
- sanitize successfully
- produce a valid 3D structure
- finish Stage 1 optimization without critical failure

### Failure rule
If Stage 1 returns `ok = 0`, **stop and diagnose before Stage 2**.

Do not continue proxy photophysics filtering on a structurally broken set.

---

## 6. Stage 2 — Semi-empirical photophysics proxy filter

### Purpose
Use a cheaper electronic-structure approximation to rank or down-select candidates.

### Typical tools
- xTB-derived screening
- sTDA / TDDFT-xTB-like proxies
- empirical Stokes-shift correction

### Critical interpretation rule
Stage 2 is a **proxy stage**, not the final truth stage.

Therefore:
- **do not use the final emission window directly as the Stage 2 hard cutoff**
- instead use a **broadened proxy window** that is wider than the final target window

### Default proxy-window rule
Unless project-specific calibration says otherwise:
- Stage-2 proxy window width should be roughly **2×** the final target-window width

### Current blue-emitter working rule
For blue TADF screening in this project:
- final target window: `450–490 nm`
- Stage-2 proxy window: `350–500 nm`
- practical oscillator-strength threshold: `f >= 0.05`

### Why this matters
Semi-empirical filters can shift relative to final DFT/TDDFT.
If Stage 2 is too strict, good candidates are discarded too early.

### Stage 2 output
A ranked shortlist with proxy photophysics fields and explicit pass/fail labels.

---

## 7. Stage 3 — Ab initio validation

### Purpose
Apply higher-level electronic-structure validation only to elite candidates.

### Typical methods
- Gaussian TDDFT
- PySCF TDDFT / state-averaged methods
- other configured ab initio engines from Stage 0

### Required protocol
For photophysical validation, use the physically meaningful order:

`S0 opt -> excited-state opt -> vertical emission`

For TADF-oriented validation, include the relevant singlet/triplet information needed for:
- `S1`
- `T1`
- `ΔE_ST` when required by the study design

### Stage 3 decision rule
Final color classification belongs here.

That means the final window such as:
- `450–490 nm` for blue

must be enforced at **Stage 3 / final validation**, not prematurely at Stage 2.

---

## 8. Stage 4 — Spectrum rendering and reporting

### Purpose
Turn validated transition information into outputs that humans can inspect or compare.

### Typical tasks
- broaden discrete transitions with Gaussian lineshapes
- convert energy to wavelength axis when needed
- export plot-ready data
- generate shortlist tables / reports / dashboard summaries

### Output examples
- spectrum tables
- emission plots
- shortlisted candidate summaries
- machine-readable CSV/JSON deliverables

---

## 9. Stage-gate policy

The workflow must obey these gates:

### Gate A — Initialization complete
If core constraints are missing, do not proceed silently.
Request clarification or emit explicit missing-field status.

### Gate B — Structural viability
If Stage 1 fails structurally, stop.
No Stage 2 photophysics filtering on invalid geometries.

### Gate C — Tool availability
If required tools are missing for a stage, emit explicit `tool_missing` status.
Do not fabricate outputs.

### Gate D — Proxy is not validation
Stage 2 may shortlist candidates, but may not claim final color confirmation.

### Gate E — Final confirmation
Only Stage 3 may claim validated photophysical outcomes.

---

## 10. Remote execution policy

When running remotely:
- record host / workdir / batch metadata
- preserve stage outputs and logs
- keep stage boundaries explicit
- do not hide partial failures behind batch wrappers

If remote and local responsibilities differ, prefer this split:
- local: library assembly, lightweight preprocessing, bookkeeping
- remote: xTB batch screening, TDDFT validation, HPC scheduling

---

## 11. Minimal canonical workflow

```text
Stage 0  Initialize workflow boundary conditions
Stage 1  Build candidates -> generate 3D -> optimize structure -> remove unstable molecules
Stage 2  Run proxy photophysics filter -> apply broadened proxy window -> create shortlist
Stage 3  Run ab initio TDDFT validation on elite candidates
Stage 4  Render spectra / export tables / summarize final conclusions
```

---

## 12. What this specification is not

This file is a **workflow contract**, not a lab notebook.

It does not try to store:
- one-off run results
- machine-specific secrets
- private infrastructure details
- ad hoc debugging logs

Those belong in run logs, reports, or operator notes.
