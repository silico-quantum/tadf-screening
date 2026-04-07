#!/usr/bin/env python3
"""
ScreeningWorkflowInitializer

A foundational initializer for generative-and-screening workflows in luminescent materials:
OLED / TADF / Fluorescence / Phosphorescence.

This version assumes NO GNN stage. Stage-1 is RDKit assembly + MMFF94 pre-optimization
followed by xTB geometry optimization/stability screening.

CRITICAL CONSTRAINT:
To estimate emission spectrum correctly, downstream steps MUST NOT use only vertical excitation from S0.
They MUST perform excited-state geometry optimization first:
- Fluorescence: S0 opt -> S1 opt -> vertical emission
- Phosphorescence: S0 opt -> T1 opt -> vertical emission
- TADF: S0 opt -> S1 opt and T1 opt (as needed) -> emission/rISC-related analysis
"""

from __future__ import annotations

import json
import math
import os
import platform
import shutil
from enum import Enum
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple

from pydantic import BaseModel, Field, model_validator


# -----------------------------
# Enums
# -----------------------------

class ComputeTier(str, Enum):
    local_basic = "local_basic"
    local_gpu = "local_gpu"
    remote_cluster = "remote_cluster"


class EmissionType(str, Enum):
    fluorescence = "Fluorescence"
    phosphorescence = "Phosphorescence"
    tadf = "TADF"


class SpectrumWidthRequirement(str, Enum):
    narrow = "Narrow/High Color Purity"
    broad = "Broad/White Light"
    any = "Any"


class TopologyType(str, Enum):
    d_a = "D-A"
    d_a_d = "D-A-D"
    a_d_a = "A-D-A"
    d_pi_a = "D-π-A"
    d_n_a = "D_n-A"


class TDDFTEngine(str, Enum):
    gaussian = "gaussian"
    orca = "orca"
    qchem = "qchem"
    pyscf = "pyscf"


# -----------------------------
# Pydantic model
# -----------------------------

class WorkflowConfig(BaseModel):
    project_name: str = Field(
        default="luminescent-screening-project",
        description="Human-readable project name for run logs and output folders.",
    )
    emission_range_nm: Tuple[float, float] = Field(
        ...,
        description="Target emission wavelength window in nm as (min_nm, max_nm), e.g. (450, 490) for blue emitters.",
    )
    emission_type: EmissionType = Field(
        ...,
        description="Emission mechanism target: Fluorescence (singlet), Phosphorescence (triplet), or TADF.",
    )
    spectrum_width_requirement: SpectrumWidthRequirement = Field(
        default=SpectrumWidthRequirement.any,
        description="Preference for spectral width/FWHM: narrow high-color-purity, broad white-light, or any.",
    )

    use_builtin_databases: List[Literal["DeepChem", "PubChem"]] = Field(
        default_factory=lambda: ["DeepChem", "PubChem"],
        description="Built-in fragment sources to enable. Supported values: DeepChem, PubChem.",
    )
    augment_with_custom_db: Optional[bool] = Field(
        default=None,
        description="Whether to augment built-in donor/acceptor fragment databases with user proprietary files.",
    )
    custom_db_paths: List[str] = Field(
        default_factory=list,
        description="Optional local paths to custom fragment libraries (.csv or .smi). Each file should contain donor/acceptor fragment SMILES records.",
    )

    topology_preferences: List[TopologyType] = Field(
        default_factory=lambda: [TopologyType.d_a, TopologyType.d_a_d],
        description="Allowed assembly topologies: D-A, D-A-D, A-D-A, D-π-A, D_n-A.",
    )

    # Stage sizing
    stage1_library_size: int = Field(
        default=10000,
        ge=100,
        description="Stage-1 candidate count for RDKit assembly + MMFF94 + xTB stability screening.",
    )
    stage2_absorption_window_nm: Tuple[float, float] = Field(
        default=(350.0, 700.0),
        description="Broad absorption window (nm) used in Stage-2 xTB-sTDDFT filtering.",
    )
    stage2_candidate_count: int = Field(
        default=1000,
        ge=10,
        description="Number of candidates promoted from Stage-1 to Stage-2 xTB-sTDDFT.",
    )
    stage3_candidate_count: int = Field(
        default=50,
        ge=1,
        description="Final elite candidate count for full TDDFT excited-state optimization validation.",
    )

    # Stokes shift strategy
    stokes_shift_strategy: Literal["fixed", "calibration", "ask_user"] = Field(
        default="ask_user",
        description="How to determine empirical Stokes shift for Stage-2 emission proxy: fixed value, calibration set, or ask user.",
    )
    empirical_stokes_shift_ev: Optional[float] = Field(
        default=None,
        description="Fixed empirical Stokes shift in eV used when strategy='fixed'. Example: 0.5 eV.",
    )
    calibration_dataset_path: Optional[str] = Field(
        default=None,
        description="Path to small calibration set (experimental or trusted references) for estimating Stokes shift when strategy='calibration'.",
    )

    # QC protocol
    tddft_level: str = Field(
        default="b3lyp/6-31+g*",
        description="Level of theory for final TDDFT validation stage (geometry + excited-state emission protocol).",
    )
    tddft_engine: Optional[TDDFTEngine] = Field(
        default=None,
        description="TDDFT engine to use (gaussian/orca/qchem/pyscf). If not set, initializer auto-discovers and selects by compute-tier default rules.",
    )

    # Runtime control
    autonomous_mode: bool = Field(
        default=False,
        description="If True, do not wait for interactive user choices; apply smart defaults automatically.",
    )
    remote_hpc_hint: Optional[str] = Field(
        default=None,
        description="Optional user-provided remote HPC hint (hostname/profile/slurm endpoint) to force remote-cluster planning.",
    )

    @model_validator(mode="after")
    def validate_ranges(self):
        mn, mx = self.emission_range_nm
        if mn <= 0 or mx <= 0 or mn >= mx:
            raise ValueError("emission_range_nm must be positive and min < max")

        a, b = self.stage2_absorption_window_nm
        if a <= 0 or b <= 0 or a >= b:
            raise ValueError("stage2_absorption_window_nm must be positive and min < max")

        if self.stokes_shift_strategy == "fixed" and self.empirical_stokes_shift_ev is None:
            raise ValueError("empirical_stokes_shift_ev is required when stokes_shift_strategy='fixed'")

        if self.stokes_shift_strategy == "calibration" and not self.calibration_dataset_path:
            raise ValueError("calibration_dataset_path is required when stokes_shift_strategy='calibration'")

        return self


# -----------------------------
# Hardware + software discovery
# -----------------------------

def _detect_hardware() -> Dict:
    hw = {
        "platform": platform.platform(),
        "os": platform.system(),
        "cpu_cores_logical": os.cpu_count() or 1,
        "ram_gb": None,
        "gpu": {"available": False, "count": 0, "names": []},
    }

    try:
        import psutil  # type: ignore

        mem = psutil.virtual_memory().total / (1024**3)
        hw["ram_gb"] = round(mem, 2)
    except Exception:
        hw["ram_gb"] = None

    try:
        import torch  # type: ignore

        if torch.cuda.is_available():
            count = torch.cuda.device_count()
            names = [torch.cuda.get_device_name(i) for i in range(count)]
            hw["gpu"] = {"available": True, "count": count, "names": names}
    except Exception:
        pass

    return hw


def _detect_hpc_slurm(remote_hpc_hint: Optional[str] = None) -> Dict:
    has_sbatch = shutil.which("sbatch") is not None
    has_squeue = shutil.which("squeue") is not None
    env_slurm = any(k.startswith("SLURM_") for k in os.environ.keys())

    return {
        "remote_hpc_hint": remote_hpc_hint,
        "slurm_detected": bool(has_sbatch or has_squeue or env_slurm or remote_hpc_hint),
        "sbatch_in_path": has_sbatch,
        "squeue_in_path": has_squeue,
        "slurm_env_present": env_slurm,
    }


def _assign_compute_tier(hw: Dict, hpc: Dict) -> ComputeTier:
    if hpc.get("slurm_detected"):
        return ComputeTier.remote_cluster
    if hw.get("gpu", {}).get("available"):
        return ComputeTier.local_gpu
    return ComputeTier.local_basic


def _discover_tddft_engines() -> Dict[str, bool]:
    # CLI binaries
    gaussian_found = any(shutil.which(cmd) for cmd in ["g16", "g09", "gaussian"])
    orca_found = shutil.which("orca") is not None
    qchem_found = shutil.which("qchem") is not None

    # Python package
    try:
        import importlib.util

        pyscf_found = importlib.util.find_spec("pyscf") is not None
    except Exception:
        pyscf_found = False

    return {
        "gaussian": bool(gaussian_found),
        "orca": bool(orca_found),
        "qchem": bool(qchem_found),
        "pyscf": bool(pyscf_found),
    }


def _choose_tddft_engine(
    requested: Optional[TDDFTEngine],
    available_map: Dict[str, bool],
    compute_tier: ComputeTier,
    autonomous_mode: bool,
) -> Tuple[Optional[str], List[str]]:
    available = [k for k, v in available_map.items() if v]
    questions: List[str] = []

    if requested is not None:
        if available_map.get(requested.value):
            return requested.value, questions
        questions.append(
            f"Requested TDDFT engine '{requested.value}' is not available. Available: {available or ['none']}."
        )

    if not available:
        questions.append("No supported TDDFT engine found (gaussian/orca/qchem/pyscf). Please install one.")
        return None, questions

    # If multiple available and not autonomous, ask user
    if len(available) > 1 and not autonomous_mode and requested is None:
        questions.append(f"Multiple TDDFT engines detected: {available}. Which one should be used?")

    # Smart defaults
    if compute_tier in (ComputeTier.local_basic, ComputeTier.local_gpu):
        if "pyscf" in available:
            return "pyscf", questions
    if compute_tier == ComputeTier.remote_cluster:
        if "gaussian" in available:
            return "gaussian", questions

    # Fallback first available
    return available[0], questions


# -----------------------------
# Fragment DB + topology helpers
# -----------------------------

def _validate_custom_db_paths(paths: List[str]) -> Dict:
    report = {"valid": [], "invalid": []}
    for p in paths:
        pp = Path(p).expanduser()
        if not pp.exists():
            report["invalid"].append({"path": p, "reason": "not_found"})
            continue
        if pp.suffix.lower() not in {".csv", ".smi"}:
            report["invalid"].append({"path": p, "reason": "unsupported_extension"})
            continue
        report["valid"].append(str(pp))
    return report


def count_halogen_leaving_groups(smiles: str) -> int:
    """Count halogens [Cl,Br,I] as proxy reactive leaving groups.
    RDKit preferred; regex fallback if RDKit unavailable.
    """
    try:
        from rdkit import Chem  # type: ignore

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0
        count = 0
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in {"Cl", "Br", "I"}:
                count += 1
        return count
    except Exception:
        # coarse fallback
        return smiles.count("Cl") + smiles.count("Br") + smiles.count("I")


def validate_topology_feasibility(topology: TopologyType, donor_smiles: str, acceptor_smiles: str) -> Tuple[bool, str]:
    """Minimal feasibility checks before RDKit assembly.

    Constraint from user:
    - D-A-D requires acceptor with >=2 leaving groups.
    """
    if topology == TopologyType.d_a_d:
        n = count_halogen_leaving_groups(acceptor_smiles)
        if n < 2:
            return False, f"skip: topology D-A-D requires acceptor with >=2 leaving groups, found {n}"
    if topology == TopologyType.a_d_a:
        n = count_halogen_leaving_groups(donor_smiles)
        if n < 2:
            return False, f"skip: topology A-D-A requires donor with >=2 leaving groups, found {n}"
    return True, "ok"


def assemble_with_rdkit_template_note() -> Dict:
    """Return a declarative note for downstream assembly implementation constraints.

    Required methods:
    - rdkit.Chem.AllChem.ReactionFromSmarts
    - rdkit.Chem.ReplaceSubstructs
    """
    return {
        "required_rdkit_assembly_methods": [
            "ReactionFromSmarts",
            "ReplaceSubstructs",
        ],
        "allowed_topologies": [t.value for t in TopologyType],
        "validation_rule": "Before D-A-D assembly, verify acceptor has >=2 halogen leaving groups [Cl,Br,I].",
    }


# -----------------------------
# Spectrum broadening
# -----------------------------

def _fwhm_to_sigma(fwhm_ev: float) -> float:
    return fwhm_ev / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def generate_emission_spectrum(
    energies_ev: List[float],
    oscillator_strengths: List[float],
    width_mode: SpectrumWidthRequirement = SpectrumWidthRequirement.any,
    fwhm_ev: Optional[float] = None,
    n_grid: int = 2000,
) -> Dict:
    """Generate broadened emission curve and return JSON-ready data.

    I(E) = sum_i f_i * exp(-(E - E_i)^2 / (2 sigma^2))
    Convert X-axis to wavelength using lambda(nm)=1240/E.
    """
    import numpy as np  # type: ignore

    if len(energies_ev) != len(oscillator_strengths):
        raise ValueError("energies_ev and oscillator_strengths must have same length")

    if not energies_ev:
        return {"energies_ev": [], "wavelength_nm": [], "intensity": [], "sticks": []}

    if fwhm_ev is None:
        if width_mode == SpectrumWidthRequirement.narrow:
            fwhm_ev = 0.15
        elif width_mode == SpectrumWidthRequirement.broad:
            fwhm_ev = 0.30
        else:
            fwhm_ev = 0.22

    sigma = _fwhm_to_sigma(fwhm_ev)

    e_min = max(0.5, min(energies_ev) - 1.0)
    e_max = max(1.0, max(energies_ev) + 1.0)
    grid_e = np.linspace(e_min, e_max, n_grid)
    intensity = np.zeros_like(grid_e)

    for e, f in zip(energies_ev, oscillator_strengths):
        intensity += f * np.exp(-((grid_e - e) ** 2) / (2.0 * sigma**2))

    # Convert axis to wavelength; keep monotonic by reversing if needed.
    wl_nm = 1240.0 / grid_e

    return {
        "meta": {
            "fwhm_ev": fwhm_ev,
            "sigma_ev": sigma,
            "width_mode": width_mode.value,
        },
        "energy_axis_ev": grid_e.tolist(),
        "wavelength_nm": wl_nm.tolist(),
        "intensity": intensity.tolist(),
        "sticks": [
            {"energy_ev": e, "wavelength_nm": 1240.0 / e, "oscillator_strength": f}
            for e, f in zip(energies_ev, oscillator_strengths)
            if e > 0
        ],
    }


# -----------------------------
# Main initializer
# -----------------------------

def initialize_workflow(params: WorkflowConfig) -> str:
    """Initialize and summarize workflow setup as a JSON string.

    Returns:
        JSON string with project boundary conditions, hardware limits,
        stage-by-stage filtering protocol, and strict downstream constraints.
    """

    hw = _detect_hardware()
    hpc = _detect_hpc_slurm(params.remote_hpc_hint)
    compute_tier = _assign_compute_tier(hw, hpc)

    db_report = _validate_custom_db_paths(params.custom_db_paths)
    engine_map = _discover_tddft_engines()
    chosen_engine, engine_questions = _choose_tddft_engine(
        requested=params.tddft_engine,
        available_map=engine_map,
        compute_tier=compute_tier,
        autonomous_mode=params.autonomous_mode,
    )

    questions: List[str] = []
    questions.extend(engine_questions)

    if params.augment_with_custom_db is None:
        questions.append(
            "Would you like to augment built-in DeepChem/PubChem donor-acceptor fragment sets with your own .csv/.smi libraries?"
        )

    if params.stokes_shift_strategy == "ask_user" and params.empirical_stokes_shift_ev is None:
        questions.append(
            "Would you like to use a fixed empirical Stokes Shift (e.g., 0.5 eV) or provide a small calibration set to estimate it?"
        )

    # Emission protocol by type (strict)
    if params.emission_type == EmissionType.fluorescence:
        excited_opt_state = "S1"
    elif params.emission_type == EmissionType.phosphorescence:
        excited_opt_state = "T1"
    else:
        excited_opt_state = "S1/T1"

    config = {
        "project": {
            "name": params.project_name,
            "compute_tier": compute_tier.value,
            "autonomous_mode": params.autonomous_mode,
        },
        "system_hardware_profile": {
            "hardware": hw,
            "hpc_slurm": hpc,
            "resource_guidance": {
                "local_basic": "Keep expensive excited-state jobs minimal; prioritize Stage-1/2 filtering.",
                "local_gpu": "Use GPU for optional acceleration, but Stage-1 remains RDKit+MMFF94+xTB stability screening.",
                "remote_cluster": "Use remote scheduler for Stage-2/3 and larger candidate throughput.",
            },
        },
        "photophysical_targets": {
            "emission_range_nm": list(params.emission_range_nm),
            "emission_type": params.emission_type.value,
            "spectrum_width_requirement": params.spectrum_width_requirement.value,
            "default_fwhm_ev": (
                0.15
                if params.spectrum_width_requirement == SpectrumWidthRequirement.narrow
                else 0.30
                if params.spectrum_width_requirement == SpectrumWidthRequirement.broad
                else 0.22
            ),
        },
        "fragment_database": {
            "builtin_sources": params.use_builtin_databases,
            "augment_with_custom_db": params.augment_with_custom_db,
            "custom_db_paths": params.custom_db_paths,
            "custom_db_validation": db_report,
            "ingestion_note": "Custom fragment files must be .csv or .smi and should include donor/acceptor SMILES records.",
        },
        "assembly_topology_rules": {
            "topology_preferences": [t.value for t in params.topology_preferences],
            "rdkit_splicing_constraints": assemble_with_rdkit_template_note(),
        },
        "quantum_chemistry_protocol": {
            "tddft_level": params.tddft_level,
            "tddft_engine_selected": chosen_engine,
            "tddft_engine_discovery": engine_map,
            "critical_emission_constraint": (
                "Do NOT estimate emission from S0-only vertical excitation. "
                f"Must optimize excited-state geometry ({excited_opt_state}) before vertical emission energy calculation."
            ),
            "strict_sequence": [
                "S0 geometry optimization",
                f"{excited_opt_state} geometry optimization",
                "Vertical emission energy calculation",
            ],
        },
        "pyramid_filtering_strategy": {
            "stage1_rdkit_mmff_xtb_stability": {
                "objective": "Assemble molecules by topology rules, run RDKit MMFF94 pre-optimization, then xTB geometry optimization/stability screening.",
                "library_size": params.stage1_library_size,
                "required_steps": [
                    "RDKit SMILES sanitization",
                    "RDKit 3D embedding",
                    "AllChem.MMFFOptimizeMolecule",
                    "xtb <file>.xyz --opt",
                ],
                "output": "Structurally valid and stable candidates promoted to Stage-2.",
            },
            "stage2_xtb_stddft_empirical_shift": {
                "objective": "xTB/sTDDFT proxy + empirical Stokes shift correction.",
                "input_count": params.stage2_candidate_count,
                "xtb_commands": [
                    "xtb <file>.xyz --opt",
                    "xtb <file>.xyz --stda",
                ],
                "stokes_shift_strategy": params.stokes_shift_strategy,
                "empirical_stokes_shift_ev": params.empirical_stokes_shift_ev,
                "calibration_dataset_path": params.calibration_dataset_path,
                "absorption_window_nm": list(params.stage2_absorption_window_nm),
                "filter_logic": "lambda_em_xtb proxy from absorption with empirical shift; retain candidates near target emission window.",
            },
            "stage3_full_tddft_validation": {
                "objective": "Final high-fidelity TDDFT validation on elite candidates.",
                "input_count": params.stage3_candidate_count,
                "protocol": "S0 opt -> excited-state opt (S1 or T1) -> vertical emission",
                "level_of_theory": params.tddft_level,
                "engine": chosen_engine,
            },
        },
        "software_stack_constraints": {
            "phase1_generation_preopt": {
                "required_package": "rdkit",
                "required_operations": [
                    "Chem.MolFromSmiles (SMILES parsing/sanitization)",
                    "AllChem.EmbedMolecule (3D conformer embedding)",
                    "AllChem.MMFFOptimizeMolecule (MMFF94 pre-optimization)",
                ],
            },
            "phase2_3_semiempirical": {
                "required_software": "xtb",
                "required_operations": [
                    "xtb <file>.xyz --opt (GFN2-xTB geometry optimization)",
                    "xtb <file>.xyz --stda (sTD-DFT excitation proxy)",
                    "Parse energies/gradients/wavelengths from text output",
                ],
            },
            "phase4_tddft_discovery": {
                "supported_engines": ["gaussian", "orca", "qchem", "pyscf"],
                "default_logic": {
                    "local_basic/local_gpu": "pyscf",
                    "remote_cluster": "gaussian",
                },
            },
            "phase5_emission_spectrum": {
                "method": "Gaussian broadening from discrete transitions",
                "equation": "I(E)=sum_i f_i exp(-(E-E_i)^2/(2 sigma^2))",
                "axis_conversion": "lambda(nm)=1240/E(eV)",
            },
        },
        "agent_questions_if_missing": questions,
    }

    return json.dumps(config, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    # Minimal demo usage
    demo = WorkflowConfig(
        emission_range_nm=(450, 490),
        emission_type=EmissionType.tadf,
        spectrum_width_requirement=SpectrumWidthRequirement.narrow,
    )
    print(initialize_workflow(demo))
