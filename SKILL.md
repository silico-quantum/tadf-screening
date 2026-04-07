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
4. Generate candidate structures (`scripts/generate_10k_structures_fast.py`).
5. Run xTB batch pre-screening (`scripts/run_xtb_batch_manifest.py`).
6. Optionally run TDDFT-xTB wavelength filter (`scripts/run_tddft_xtb_filter.py`).
7. Send elite candidates to full TDDFT validation.

Do not skip stage gates.

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

## Remote execution marking (mandatory)

Before starting each screening run, update `BATCH_STATUS.md` with:
- batch id
- stage
- execution target (`local` or `remote`)
- remote host/workdir (required for remote)
- start time
- process id

## References

- Workflow spec: `references/workflow-spec.md`
- Stage I/O contract: `references/io-contract.md`
- Topology assembly script: `scripts/build_da_topology_library.py`
