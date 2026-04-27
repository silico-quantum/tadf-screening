#!/usr/bin/env python3
"""
Stage 4: MOMAP Photophysics Screening

Integrates momap-tadf into the TADF screening pipeline.
Takes Stage 3 Gaussian outputs (S0/T1/S1 logs) and runs:
  1. EVC (electron-vibration coupling) S1→S0
  2. spec_tvcf (fluorescence spectrum)
  3. ISC (S1→T1 intersystem crossing rate)
  4. Ranking by target window + oscillator strength + ΔE_ST

Usage:
  python stage4_momap.py candidates.csv --target blue
  python stage4_momap.py --mol-id mol_07566 --s0 s0.log --s1 s1.log --t1 t1.log --target green

Color presets (--target):
  deep-blue  430-460 nm    sky-blue  470-500 nm
  blue       450-490 nm    green     500-550 nm
  yellow     550-600 nm    red       600-700 nm
  nir        700-1000 nm   custom    (use --window MIN MAX)
"""
import sys
import os
import csv
import json
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

# Color emission window presets (nm)
COLOR_PRESETS = {
    'deep-blue': (430, 460),
    'blue':      (450, 490),
    'sky-blue':  (470, 500),
    'green':     (500, 550),
    'yellow':    (550, 600),
    'orange':    (580, 620),
    'red':       (600, 700),
    'deep-red':  (650, 750),
    'nir':       (700, 1000),
}

# Path to momap tools
TOOLS_DIR = os.path.join(os.path.dirname(__file__), '..', 'momap', 'tools')
if not os.path.isdir(TOOLS_DIR):
    # Try absolute fallback
    TOOLS_DIR = os.path.expanduser('~/.openclaw/skills/momap/tools')

def run_single(mol_id, s0_log, s1_log, t1_log, output_dir, temperature=300):
    """Run MOMAP TADF pipeline for one molecule via tadf.py."""
    tadf_script = os.path.join(TOOLS_DIR, 'tadf.py')
    if not os.path.isfile(tadf_script):
        # Import directly
        sys.path.insert(0, TOOLS_DIR)
        from tadf import process_molecule
        return process_molecule(mol_id, s0_log, s1_log, t1_log, output_dir, temperature)

    cmd = [
        sys.executable, tadf_script, mol_id,
        '--s0', s0_log,
        '--s1', s1_log,
        '--output', output_dir,
        '--temperature', str(temperature),
        '--json',
    ]
    if t1_log and Path(t1_log).exists():
        cmd.extend(['--t1', t1_log])

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if result.returncode == 0:
        return json.loads(result.stdout)
    else:
        print(f"  ❌ {mol_id}: {result.stderr[-200:]}")
        return {'mol_id': mol_id, 'success': False, 'error': result.stderr[-200:]}

def rank_candidates(results, target_window=(450, 490)):
    """Rank TADF candidates by target window + ΔE_ST + oscillator strength."""
    win_min, win_max = target_window
    for r in results:
        if not r.get('success'):
            continue
        
        score = 0
        
        # Target window bonus
        spec = r.get('spectrum', {})
        peak_nm = spec.get('peak_wavelength', 0)
        if win_min <= peak_nm <= win_max:
            score += 3.0  # inside target window
        elif win_min - 10 <= peak_nm <= win_max + 10:
            score += 1.5  # near window (10 nm margin)
        elif win_min - 30 <= peak_nm <= win_max + 30:
            score += 0.5  # visible vicinity
        
        # ΔE_ST (smaller is better for TADF)
        delta_est = r.get('delta_EST_eV', 1.0)
        if delta_est and delta_est > 0:
            if delta_est < 0.1:
                score += 2.0
            elif delta_est < 0.2:
                score += 1.5
            elif delta_est < 0.3:
                score += 1.0
            elif delta_est < 0.5:
                score += 0.5
        
        # Oscillator strength (higher is better)
        f_emi = r.get('f_emi', 0)
        if f_emi > 0.1:
            score += 1.0
        elif f_emi > 0.05:
            score += 0.7
        elif f_emi > 0.01:
            score += 0.3
        
        r['_score'] = round(score, 2)
        r['_in_window'] = (win_min <= peak_nm <= win_max)
        r['_near_window'] = (win_min - 10 <= peak_nm <= win_max + 10)
    
    return sorted(results, key=lambda r: r.get('_score', -1), reverse=True)


def generate_report(results, output_path, target_window=(450, 490), target_name='blue'):
    """Generate a markdown photophysics report."""
    ranked = rank_candidates(results, target_window)
    win_min, win_max = target_window
    
    lines = [
        f"# TADF Stage 4 — MOMAP Photophysics Report",
        f"",
        f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"**Total candidates:** {len(results)}",
        f"**Successfully processed:** {sum(1 for r in results if r.get('success'))}",
        f"",
        f"## 🏆 Top Candidates (Ranked)",
        f"",
        f"| Rank | Molecule | λ_emi (nm) | ΔE_ST (eV) | f_emi | Ead (eV) | {target_name.capitalize()}? | Score |",
        f"|------|----------|-----------|-----------|-------|---------|-------|-------|",
    ]
    
    for i, r in enumerate(ranked):
        if not r.get('success'):
            continue
        if i >= 20:
            break
        
        spec = r.get('spectrum', {})
        peak = spec.get('peak_wavelength', '—')
        f_peak = f"{peak:.0f}" if isinstance(peak, (int, float)) else str(peak)
        
        delta_est = r.get('delta_EST_eV')
        dE_str = f"{delta_est:.3f}" if delta_est else "—"
        
        f_emi = r.get('f_emi', 0)
        f_str = f"{f_emi:.4f}" if f_emi else "—"
        
        ead = r.get('Ead_S1_S0')
        ead_str = f"{ead:.3f}" if ead else "—"
        
        bw = r.get('_in_window')
        near = r.get('_near_window')
        if bw:
            blue_str = '✅ IN'
        elif near:
            blue_str = '🟡 near'
        else:
            blue_str = '⚪'
        
        score = r.get('_score', 0)
        
        lines.append(f"| {i+1} | {r['mol_id']} | {f_peak} | {dE_str} | {f_str} | {ead_str} | {blue_str} | {score:.1f} |")
    
    lines.extend([
        "",
        "## 📊 Screening Summary",
        "",
    ])
    
    blue_count = sum(1 for r in ranked if r.get('_in_window'))
    near_count = sum(1 for r in ranked if r.get('_near_window') and not r.get('_in_window'))
    success_count = sum(1 for r in ranked if r.get('success'))
    
    lines.append(f"- ✅ **{target_name} window ({win_min}–{win_max} nm):** {blue_count}/{success_count}")
    lines.append(f"- 🟡 **Near {target_name} (±10 nm):** {near_count}/{success_count}")
    lines.append(f"- ✅ **Successfully processed:** {success_count}/{len(results)}")
    
    if blue_count > 0:
        lines.append(f"\n### ✅ {target_name.capitalize()} Window Candidates")
        for r in ranked:
            if r.get('_in_window'):
                spec = r.get('spectrum', {})
                peak = spec.get('peak_wavelength', '—')
                lines.append(f"- **{r['mol_id']}**: λ_emi={peak:.0f} nm, ΔE_ST={r.get('delta_EST_eV', 0):.3f} eV")
    
    lines.extend([
        "",
        "---",
        f"*Stage 4 powered by MOMAP 2024A (TVCF method, B3LYP/6-31G*) | Silico 🔮*"
    ])
    
    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))
    
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description='Stage 4: MOMAP Photophysics — TADF screening finale'
    )
    parser.add_argument('input', nargs='?', help='CSV with columns: mol_id,s0_log,s1_log,t1_log')
    parser.add_argument('--mol-id', help='Single molecule ID')
    parser.add_argument('--s0', help='S0 Gaussian log (single mode)')
    parser.add_argument('--s1', help='S1 Gaussian log (single mode)')
    parser.add_argument('--t1', help='T1 Gaussian log (single mode)')
    parser.add_argument('--output', '-o', default='./stage4_output', help='Output directory')
    parser.add_argument('--temperature', '-T', type=float, default=300, help='Temperature (K)')
    parser.add_argument('--target', default='blue',
                       choices=list(COLOR_PRESETS.keys()),
                       help='Target emission color (preset window). Can also be set via --config from Stage 0.')
    parser.add_argument('--window', '-w', nargs=2, type=float, metavar=('MIN', 'MAX'),
                       help='Custom emission window in nm (e.g. --window 480 520)')
    parser.add_argument('--config', '-c', help='Path to Stage 0 workflow config JSON (auto-extracts target color/window)')
    parser.add_argument('--report', '-r', default='stage4_report.md', help='Report filename')
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Resolve target window
    if args.config and Path(args.config).exists():
        with open(args.config) as f:
            cfg = json.load(f)
        targets = cfg.get('photophysical_targets', {})
        cfg_color = targets.get('emission_color', '')
        cfg_range = targets.get('emission_range_nm', [])
        if cfg_color and cfg_color in COLOR_PRESETS:
            target_window = COLOR_PRESETS[cfg_color]
            target_name = cfg_color
        elif len(cfg_range) == 2:
            target_window = tuple(cfg_range)
            target_name = f'{target_window[0]}-{target_window[1]}'
        else:
            target_window = COLOR_PRESETS['blue']
            target_name = 'blue'
        print(f"📋 Loaded from config: {target_name} ({target_window[0]}–{target_window[1]} nm)")
    elif args.window:
        target_window = tuple(args.window)
        target_name = f'{target_window[0]}-{target_window[1]}'
    else:
        target_window = COLOR_PRESETS.get(args.target, COLOR_PRESETS['blue'])
        target_name = args.target
    
    print(f"🎯 Target: {target_name} ({target_window[0]}–{target_window[1]} nm)")

    results = []

    # Single molecule mode
    if args.mol_id:
        if not args.s0 or not args.s1:
            print("❌ Single mode requires --s0 and --s1")
            return 1
        
        result = run_single(args.mol_id, args.s0, args.s1, args.t1, str(output_dir), args.temperature)
        results.append(result)
    
    # Batch CSV mode
    elif args.input:
        with open(args.input) as f:
            reader = csv.DictReader(f)
            for row in reader:
                mol_id = row.get('mol_id', '').strip()
                if not mol_id:
                    continue
                
                s0 = row.get('s0_log', '').strip()
                s1 = row.get('s1_log', '').strip()
                t1 = row.get('t1_log', '').strip() if 't1_log' in row else None
                
                if not s0 or not s1:
                    print(f"  ⚠️  {mol_id}: missing S0 or S1 log, skipping")
                    continue
                
                result = run_single(mol_id, s0, s1, t1, str(output_dir), args.temperature)
                results.append(result)
    
    else:
        parser.print_help()
        return 1

    # Generate report
    if results:
        report_path = output_dir / args.report
        generate_report(results, report_path, target_window, target_name)
        print(f"\n📋 Report: {report_path}")
        
        # Also dump JSON
        json_path = output_dir / 'stage4_results.json'
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"📊 JSON: {json_path}")

        # Quick summary
        ranked = rank_candidates(results, target_window)
        in_window = [r for r in ranked if r.get('_in_window')]
        print(f"\n{'='*50}")
        print(f"🏁 Stage 4 Complete")
        print(f"{'='*50}")
        print(f"  Candidates: {len(results)}")
        print(f"  Successful: {sum(1 for r in results if r.get('success'))}")
        print(f"  In {target_name} window: {len(in_window)}")
        if in_window:
            for r in in_window[:5]:
                spec = r.get('spectrum', {})
                print(f"    🔵 {r['mol_id']}: {spec.get('peak_wavelength', '?'):.0f} nm  ΔE_ST={r.get('delta_EST_eV', 0):.3f} eV  score={r.get('_score', 0):.1f}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
