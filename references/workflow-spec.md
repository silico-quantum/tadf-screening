# Workflow Specification

## Stage 0: Initialization
Run `scripts/screening_workflow_initializer.py` to determine:
- compute tier
- target emission window/type
- database boundaries
- TDDFT engine and level of theory
- inquiry-stage runtime mode

## Stage 1: Structure generation and stability pre-screen
- Assemble molecules with topology rules.
- Generate 3D geometry with RDKit.
- Run MMFF94 pre-optimization.
- Export `.xyz`.
- Run `xtb --opt` for structural stability.

## Stage 2: Semi-empirical photophysics filter
- Run `xtb --stda` for surviving structures.
- Apply empirical Stokes shift.
- Keep candidates inside target emission window.
- For TADF, enforce practical S1-T1 gap gate as configured.

## Stage 3: Ab initio validation
- Run full TDDFT validation on elite candidates.
- Mandatory order: `S0 opt -> excited-state opt -> vertical emission`.

## Stage 4: Emission spectrum rendering
- Use Gaussian broadening on discrete transitions.
- Convert energy axis to wavelength axis.
- Export plot-ready or dashboard-ready data.
