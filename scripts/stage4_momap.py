#!/usr/bin/env python3
"""
Stage 4: MOMAP Photophysics Screening

Integrates momap-tadf into the TADF screening pipeline.
Takes Stage 3 Gaussian outputs (S0/T1/S1 logs) and runs:
  1. EVC (electron-vibration coupling) S1→S0
  2. spec_tvcf (fluorescence spectrum)
  3. ISC (S1→T1 intersystem crossing rate)
  4. Ranking by blue window (450–490 nm) + oscillator strength + ΔE_ST

Usage:
  python stage4_momap.py candidates.csv --output results/
  python stage4_momap.py --mol-id mol_07566 --s0 s0.log --s1 s1.log --t1 t1.log
"""
import sys
import os
import csv
import json
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

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

def rank_candidates(results):
    """Rank TADF candidates by blue window + ΔE_ST + oscillator strength."""
    for r in results:
        if not r.get('success'):
            continue
        
        score = 0
        
        # Blue window bonus
        spec = r.get('spectrum', {})
        peak_nm = spec.get('peak_wavelength', 0)
        if 450 <= peak_nm <= 490:
            score += 3.0  # inside blue window
        elif 440 <= peak_nm <= 500:
            score += 1.5  # near blue
        elif 400 <= peak_nm <= 550:
            score += 0.5  # visible
        
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
    
    return sorted(results, key=lambda r: r.get('_score', -1), reverse=True)


def generate_report(results, output_path):
    """Generate a markdown photophysics report."""
    ranked = rank_candidates(results)
    
    lines = [
        f"# TADF Stage 4 — MOMAP Photophysics Report",
        f"",
        f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"**Total candidates:** {len(results)}",
        f"**Successfully processed:** {sum(1 for r in results if r.get('success'))}",
        f"",
        f"## 🏆 Top Candidates (Ranked)",
        f"",
        f"| Rank | Molecule | λ_emi (nm) | ΔE_ST (eV) | f_emi | Ead (eV) | Blue? | Score |",
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
        
        bw = r.get('blue_window')
        if bw == True:
            blue_str = '🔵 YES'
        elif bw == 'near':
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
    
    blue_count = sum(1 for r in ranked if r.get('blue_window') == True)
    near_count = sum(1 for r in ranked if r.get('blue_window') == 'near')
    success_count = sum(1 for r in ranked if r.get('success'))
    
    lines.append(f"- 🔵 **Blue window (450–490 nm):** {blue_count}/{success_count}")
    lines.append(f"- 🟡 **Near blue (440–500 nm):** {near_count}/{success_count}")
    lines.append(f"- ✅ **Successfully processed:** {success_count}/{len(results)}")
    
    if blue_count > 0:
        lines.append(f"\n### 🔵 Blue Window Candidates")
        for r in ranked:
            if r.get('blue_window') == True:
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
    parser.add_argument('--report', '-r', default='stage4_report.md', help='Report filename')
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

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
        generate_report(results, report_path)
        print(f"\n📋 Report: {report_path}")
        
        # Also dump JSON
        json_path = output_dir / 'stage4_results.json'
        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"📊 JSON: {json_path}")

        # Quick summary
        ranked = rank_candidates(results)
        blue = [r for r in ranked if r.get('blue_window') == True]
        print(f"\n{'='*50}")
        print(f"🏁 Stage 4 Complete")
        print(f"{'='*50}")
        print(f"  Candidates: {len(results)}")
        print(f"  Successful: {sum(1 for r in results if r.get('success'))}")
        print(f"  Blue window: {len(blue)}")
        if blue:
            for r in blue[:5]:
                spec = r.get('spectrum', {})
                print(f"    🔵 {r['mol_id']}: {spec.get('peak_wavelength', '?'):.0f} nm  ΔE_ST={r.get('delta_EST_eV', 0):.3f} eV  score={r.get('_score', 0):.1f}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
