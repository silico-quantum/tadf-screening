---
name: tadf-screening
description: Standardized high-throughput screening workflow for luminescent materials (TADF, Fluorescence, Phosphorescence). Use when the user asks to generate candidates, run xTB pre-screening, apply TDDFT-xTB wavelength filtering, or launch final TDDFT validation on Marcus/remote clusters.
---

# TADF Screening Skill (Standard)

Run this skill for end-to-end candidate generation and filtering.

## Required execution order

1. Run `scripts/screening_workflow_initializer.py` first.
2. Confirm topology and initial sample count during inquiry stage (default initial sample count: `10000`).
3. Build topology-driven D/A assemblies (`scripts/build_da_topology_library.py`).
4. Generate candidate 3D structures with your project-level generator (outside this skill package) and export `.xyz` files.
5. Run xTB batch pre-screening (`scripts/run_xtb_batch_manifest.py`).
6. Optionally run TDDFT-xTB wavelength filter (`scripts/run_tddft_xtb_filter.py`).
7. Send elite candidates to full TDDFT validation.

Do not skip stage gates.

## Supervision requirement (all stages)

Every stage must include an explicit supervision task and report item.

- **Stage 0 (initializer)**: verify critical inquiry fields are resolved; report unresolved blockers.
- **Stage 1 (SMILES -> XYZ)**: track `xyz_count / target_count`, manifest existence, process alive/stopped, restart action if interrupted.
- **Stage 2 (xTB + sTDA)**: track Slurm queue state, produced `_opt.log` / `_xtb4stda.log` / `_stda.log` counts, and error signature checks (e.g., invalid option errors).
- **Stage 3 (TDDFT validation)**: track submitted jobs, running/pending/failed counts, parsed success ratio, and failed-case reasons.
- **Stage 4 (analysis/report)**: track result file completeness (ranked CSV, shortlist, plots) and publish final summary with pass/fail statistics.

Minimum monitoring cadence default remains:
- check every 5 minutes
- report every 10 minutes

Timed report content (mandatory):
- **Compute process status**: running/pending/completed/failed counts (or alive/stopped for local stages), plus scheduler job IDs when remote.
- **File generation status**: key output file counts (e.g., logs, result CSVs, optimized structures) versus expected totals.
- **Error status**: whether errors are present; include latest error signature/snippet and recovery action if any.

When monitoring is active, every periodic report must include all three items above in compact form.

## 3D structure generation: multiprocessing rule

When the initial sample count exceeds **10,000**, the 3D generation stage MUST use multiprocessing:

- **Auto-enable**: `--workers -1` (default) → `cpu_count` workers when target > 10k
- **Override**: `--workers N` to set explicit worker count; `--workers 0` for single-process
- **Assembly is sequential** (fast); only the heavy ETKDG + MMFF94 step is parallelized
- Batch assembly size: 50,000 SMILES, then dispatch to worker pool
- Supports resume via checkpoint.json

Rationale: single-process ETKDG yields ~1,200 mol/min on ARM64; multiprocessing with 8-10 workers can reach ~8,000-10,000 mol/min.

## Hard chemistry constraint

For emission prediction, do **not** rely on S0-only vertical excitation.
Downstream validation must use:
- `S0 optimization -> excited-state optimization (S1 for Fluorescence/TADF or T1 for Phosphorescence) -> vertical emission calculation`.

## Inquiry-stage runtime mode (must be set)

The initializer supports inquiry-stage operation control:
- `interaction_mode`: `interactive_chat | config_only | autonomous`
- `question_mode`: `blocking | non_blocking`
- `max_question_rounds`: integer >= 1

Policy:
- In `autonomous` mode, use `non_blocking` question mode.
- Record all unresolved questions in initializer output under `agent_questions_if_missing`.
- Active inquiry is mandatory: if critical fields are not provided (emission range/type, spectrum width preference, empirical Stokes shift, topology), ask explicitly before proceeding in `blocking` mode.
- Monitoring cadence is also inquiry-driven: ask user for check/report intervals. Defaults are check every 5 minutes and report every 10 minutes.
- **Resource inquiry must run first**: confirm remote resource profile, check `xtb4stda/stda` availability, and explicitly ask user permission before attempting installation of missing tools.

## Topology and assembly constraints

Use only rule-based assembly templates:
- `D-A`, `D-A-D`, `A-D-A`, `D-pi-A`, `D_n-A`
- Use RDKit `ReactionFromSmarts` or `ReplaceSubstructs`
- Before `D-A-D`, verify acceptor has >=2 leaving groups (`Cl`, `Br`, `I`)

Use `scripts/build_da_topology_library.py` to enforce topology-specific assembly.
The inquiry stage must explicitly confirm:
- selected topology (single or mixed)
- initial sample count (default `10000`)

## xTB and filtering contract

### sTDA invocation standard (critical)

- Do **not** call `xtb --stda` (unsupported in current xTB build and will fail).
- Use the two-step flow:
  1. `xtb4stda <optimized.xyz>`
  2. `stda -xtb -e <N>` (example: `stda -xtb -e 10`)
- If `xtb4stda` or `stda` is missing, stop with explicit `tool_missing` and do not continue silently.

Input manifest must include:
```csv
idx,name,xyz_path
```

xTB output contract (`xtb_progress.csv`) must include:
```csv
idx,name,xyz_path,status,detail,total_energy_eh,homo_lumo_gap_ev,normal_termination
```

Gate rule:
- If xTB `ok=0`, stop and diagnose before advancing.

## Local resource profiles standard (mandatory)

Use `local/` to store per-user environment settings.

- folder: `skills/tadf-screening/local/`
- template: `local/resource_profile.template.yaml`
- guidance: `local/README.md`

Rules:
- Never hardcode user-specific hosts/modules/paths in core scripts.
- Select profile during inquiry stage before remote execution.
- Record selected profile in `BATCH_STATUS.md`.

## Cron-based monitoring standard (mandatory)

Use cron for autonomous monitoring and reporting during long runs.

Default policy:
- process check interval: 5 minutes
- user report interval: 10 minutes

Use script:
- `scripts/configure_progress_cron.py --check-min 5 --report-min 10 --replace`

If user specifies other intervals, use user values.

## Local vs Remote execution boundary (mandatory)

Execution boundary for this project:
- **Local-only stage**: molecular structure file generation (`SMILES -> .xyz`, Stage-1 generation step).
- **Remote stages**: all downstream screening after `.xyz` exists (xTB batch filtering, TDDFT-xTB filter, final TDDFT validation) should run on a user-confirmed remote environment.

Resource selection policy:
- Always ask the user which remote resource/profile to use before remote execution.
- Do not hardcode a specific cluster/host as a universal default in this skill.
- Resolve host/path/module details from `local/*.yaml` after user confirmation.

## Remote execution marking (mandatory)

Before starting each screening run, update `BATCH_STATUS.md` with:
- batch id
- stage
- execution target (`local` or `remote`)
- remote host/workdir (required for remote)
- start time
- process id
- selected local resource profile (when available)

## References

- Workflow spec: `references/workflow-spec.md`
- Stage I/O contract: `references/io-contract.md`
- Topology assembly script: `scripts/build_da_topology_library.py`
- Cron setup script: `scripts/configure_progress_cron.py`
