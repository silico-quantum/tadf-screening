#!/usr/bin/env python3
"""
momap-plot — Spectrum visualization from MOMAP spec.tvcf.spec.dat

Generates publication-quality PNG with dual panels (energy + wavelength domain).
Blue window (450–490 nm) highlighting for TADF screening.
"""
import sys
import os
import argparse
from pathlib import Path

def load_spec_data(filepath):
    """Load MOMAP spec.tvcf.spec.dat."""
    data = []
    with open(filepath) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 7:
                data.append([float(x) for x in parts[:7]])
    if not data:
        raise ValueError(f"No data in {filepath}")
    
    return {
        'energy_ha':   [r[0] for r in data],
        'energy_ev':   [r[1] for r in data],
        'wavenumber':  [r[2] for r in data],
        'wavelength':  [r[3] for r in data],
        'fc_abs':      [r[4] for r in data],
        'fc_emi':      [r[5] for r in data],
        'fc_emi_int':  [r[6] for r in data],
    }

def normalize_positive(arr):
    pos = [max(v, 0) for v in arr]
    mx = max(pos) if max(pos) > 0 else 1
    return [v / mx for v in pos], mx

def plot_with_pillow(spec, output_path, title='MOMAP Spectrum', blue_window=(450, 490)):
    """Render spectrum using Pillow (no matplotlib required)."""
    from PIL import Image, ImageDraw
    
    W, H = 1600, 1000
    img = Image.new('RGB', (W, H), (12, 12, 30))
    draw = ImageDraw.Draw(img)
    
    # Colors
    EMISSION_COLOR = (79, 195, 247)   # cyan-blue
    ABSORPTION_COLOR = (255, 138, 128) # coral
    BLUE_WINDOW_COLOR = (100, 180, 255, 60)
    GRID_COLOR = (35, 35, 55)
    TEXT_COLOR = (160, 170, 190)
    TITLE_COLOR = (120, 200, 255)
    
    # Title
    draw.text((W//2 - 300, 20), title, fill=TITLE_COLOR)
    
    # === Layout ===
    margin_l, margin_r, margin_t, margin_b = 120, 60, 60, 80
    pw, ph = W - margin_l - margin_r, H - margin_t - margin_b
    half_w = pw // 2 - 20
    
    # ===== LEFT PANEL: Energy Domain =====
    lx, ly, lw, lh = margin_l, margin_t, half_w, ph
    
    ev = spec['energy_ev']
    emi_norm, emi_max = normalize_positive(spec['fc_emi'])
    abs_norm, abs_max = normalize_positive(spec['fc_abs'])
    
    # Frame & grid
    draw.rectangle([lx, ly, lx+lw, ly+lh], outline=(60, 60, 90))
    for i in range(1, 5):
        y = ly + lh - int(lh * i / 5)
        draw.line([lx, y, lx+lw, y], fill=GRID_COLOR)
        draw.text((lx-55, y-8), f"{i/5:.1f}", fill=TEXT_COLOR)
    
    draw.text((lx + lw//2 - 50, ly+lh+10), 'Energy (eV)', fill=TEXT_COLOR)
    
    # Plot emission
    ev_mask = [i for i, e in enumerate(ev) if 0.02 < e < 0.30]
    if ev_mask:
        pts = []
        for i in ev_mask:
            x = lx + int(lw * (ev[i] - 0.02) / 0.28)
            y = ly + lh - int(lh * emi_norm[i] * 0.95)
            pts.append((x, y))
        for j in range(len(pts)-1):
            draw.line([pts[j], pts[j+1]], fill=EMISSION_COLOR, width=2)
        # Fill
        for j in range(len(pts)):
            draw.line([pts[j], (pts[j][0], ly+lh)], fill=(79, 195, 247, 20))
    
    # Plot absorption
    if abs_max > 0:
        for j in range(len(pts)-1):
            y1 = ly + lh - int(lh * abs_norm[ev_mask[j]] * 0.95)
            y2 = ly + lh - int(lh * abs_norm[ev_mask[j+1]] * 0.95)
            if y1 > ly and y2 > ly:
                draw.line([(pts[j][0], y1), (pts[j+1][0], y2)], fill=ABSORPTION_COLOR, width=1)
    
    # Legend
    draw.rectangle([lx+10, ly+10, lx+180, ly+55], fill=(12, 12, 30), outline=(60, 60, 90))
    draw.line([lx+20, ly+25, lx+50, ly+25], fill=EMISSION_COLOR, width=3)
    draw.text((lx+55, ly+17), 'Emission (FC)', fill=TEXT_COLOR)
    draw.line([lx+20, ly+40, lx+50, ly+40], fill=ABSORPTION_COLOR, width=1)
    draw.text((lx+55, ly+32), 'Absorption (FC)', fill=TEXT_COLOR)
    
    # ===== RIGHT PANEL: Wavelength Domain =====
    rx, ry, rw, rh = lx + lw + 40, margin_t, half_w, ph
    
    wl = spec['wavelength']
    wl_mask = [i for i, w in enumerate(wl) if 250 < w < 800 and emi_norm[i] > 0.001]
    
    draw.rectangle([rx, ry, rx+rw, ry+rh], outline=(60, 60, 90))
    for i in range(1, 5):
        y = ry + rh - int(rh * i / 5)
        draw.line([rx, y, rx+rw, y], fill=GRID_COLOR)
    for nm in [300, 400, 500, 600, 700]:
        x = rx + int(rw * (800 - nm) / 550)
        draw.line([x, ry, x, ry+rh], fill=GRID_COLOR)
        draw.text((x-15, ry+rh+5), str(nm), fill=TEXT_COLOR)
    
    draw.text((rx + rw//2 - 50, ry+rh+10), 'Wavelength (nm)', fill=TEXT_COLOR)
    
    # Blue window highlight
    if blue_window:
        bw_start = blue_window[0]
        bw_end = blue_window[1]
        bx1 = rx + int(rw * (800 - bw_end) / 550)
        bx2 = rx + int(rw * (800 - bw_start) / 550)
        # Semi-transparent blue band
        for bx in range(bx1, bx2, 2):
            draw.line([bx, ry, bx, ry+rh], fill=(30, 60, 100, 30))
        draw.text((bx1 + 5, ry + 5), f'{bw_start}–{bw_end} nm', fill=(100, 180, 255))
    
    # Plot emission (wavelength)
    if wl_mask:
        pts = []
        for i in wl_mask:
            x = rx + int(rw * (800 - wl[i]) / 550)
            y = ry + rh - int(rh * emi_norm[i] * 0.95)
            pts.append((x, y))
        for j in range(len(pts)-1):
            draw.line([pts[j], pts[j+1]], fill=(129, 199, 132), width=2)
    
    # Peak annotations
    
    # Find top peaks
    peaks = []
    for i in range(1, len(emi_norm)-1):
        if emi_norm[i] > emi_norm[i-1] and emi_norm[i] > emi_norm[i+1] and emi_norm[i] > 0.1:
            peaks.append((wl[i], emi_norm[i]))
    peaks.sort(key=lambda x: x[1], reverse=True)
    
    for p in peaks[:3]:
        px = rx + int(rw * (800 - p[0]) / 550)
        py = ry + rh - int(rh * p[1] * 0.95)
        draw.text((px - 20, py - 20), f'{p[0]:.0f} nm', fill=(255, 220, 100))
    
    # Footer
    draw.rectangle([0, H-25, W, H], fill=(8, 8, 20))
    draw.text((W//2-280, H-22), 
              f"MOMAP 2024A | TVCF | T=300K | Peak: {peaks[0][0]:.0f} nm" if peaks else "MOMAP 2024A | TVCF | T=300K",
              fill=(100, 100, 120))
    
    img.save(output_path, 'PNG')
    return output_path

def plot_with_matplotlib(spec, output_path, title='MOMAP Spectrum', blue_window=(450, 490)):
    """Render spectrum using matplotlib (higher quality)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    
    ev = np.array(spec['energy_ev'])
    wl = np.array(spec['wavelength'])
    emi = np.maximum(np.array(spec['fc_emi']), 0)
    abs_ = np.maximum(np.array(spec['fc_abs']), 0)
    
    emi = emi / emi.max()
    if abs_.max() > 0:
        abs_ = abs_ / abs_.max()
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.patch.set_facecolor('#0c0c1e')
    
    # Energy domain
    ax = axes[0]
    ax.set_facecolor('#111133')
    mask = (ev > 0.02) & (ev < 0.30)
    ax.plot(ev[mask], emi[mask], color='#4fc3f7', lw=1.5, label='Emission')
    if abs_.max() > 0:
        ax.plot(ev[mask], abs_[mask], color='#ff8a80', lw=1, alpha=0.7, label='Absorption')
    ax.set_xlabel('Energy (eV)', color='#aaa')
    ax.set_ylabel('Normalized Intensity', color='#aaa')
    ax.legend(facecolor='#111133', edgecolor='#333', labelcolor='#ccc')
    ax.tick_params(colors='#aaa')
    ax.grid(True, alpha=0.2)
    ax.set_title('Energy Domain', color='#88c8ff')
    
    # Wavelength domain
    ax = axes[1]
    ax.set_facecolor('#111133')
    mask = (wl > 250) & (wl < 800) & (emi > 0.001)
    
    # Blue window
    bw_start, bw_end = blue_window
    ax.axvspan(bw_start, bw_end, alpha=0.15, color='#3366aa', label=f'{bw_start}–{bw_end} nm')
    
    ax.plot(wl[mask], emi[mask], color='#81c784', lw=1.5, label='Emission')
    ax.set_xlabel('Wavelength (nm)', color='#aaa')
    ax.set_ylabel('Normalized Intensity', color='#aaa')
    ax.legend(facecolor='#111133', edgecolor='#333', labelcolor='#ccc')
    ax.tick_params(colors='#aaa')
    ax.grid(True, alpha=0.2)
    ax.invert_xaxis()
    ax.set_title('Wavelength Domain', color='#81c784')
    
    # Peak annotation
    peak_idx = np.argmax(emi[mask])
    peak_nm = wl[mask][peak_idx]
    ax.annotate(f'{peak_nm:.0f} nm', xy=(peak_nm, emi[mask][peak_idx]),
               xytext=(peak_nm-30, emi[mask][peak_idx]+0.1),
               arrowprops=dict(arrowstyle='->', color='#ffd54f'),
               color='#ffd54f', fontsize=11)
    
    plt.suptitle(title, color='#ccc', fontsize=14, y=0.98)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='#0c0c1e')
    plt.close()
    return output_path

def main():
    parser = argparse.ArgumentParser(description='Plot MOMAP spectrum')
    parser.add_argument('spec_file', help='Path to spec.tvcf.spec.dat')
    parser.add_argument('--output', '-o', default=None, help='Output PNG path (default: <spec_file>.png)')
    parser.add_argument('--title', '-t', default='MOMAP TVCF Spectrum', help='Plot title')
    parser.add_argument('--blue', '-b', nargs=2, type=float, default=[450, 490],
                       help='Blue window range in nm (default: 450 490)')
    parser.add_argument('--desktop', '-D', action='store_true',
                       help='Save to Desktop')
    args = parser.parse_args()
    
    if args.output is None:
        args.output = str(Path(args.spec_file).with_suffix('.png'))
    
    if args.desktop:
        desktop = os.path.expanduser('~/Desktop')
        args.output = os.path.join(desktop, os.path.basename(args.output))
    
    spec = load_spec_data(args.spec_file)
    
    # Try matplotlib first, fall back to Pillow
    try:
        path = plot_with_matplotlib(spec, args.output, args.title, tuple(args.blue))
        print(f"✅ Spectrum saved: {path} (matplotlib)")
    except ImportError:
        path = plot_with_pillow(spec, args.output, args.title, tuple(args.blue))
        print(f"✅ Spectrum saved: {path} (Pillow)")
    
    # Print peak info
    ev = spec['energy_ev']
    wl = spec['wavelength']
    emi = [max(v, 0) for v in spec['fc_emi']]
    emi_norm = [v / max(emi) for v in emi]
    
    peaks = []
    for i in range(1, len(emi_norm)-1):
        if emi_norm[i] > emi_norm[i-1] and emi_norm[i] > emi_norm[i+1] and emi_norm[i] > 0.05:
            peaks.append({'nm': wl[i], 'ev': ev[i], 'intensity': emi_norm[i]})
    peaks.sort(key=lambda x: x['intensity'], reverse=True)
    
    print(f"\n📊 Peak Analysis:")
    for i, p in enumerate(peaks[:5]):
        marker = ''
        if args.blue and args.blue[0] <= p['nm'] <= args.blue[1]:
            marker = ' 🔵 BLUE WINDOW!'
        elif args.blue and args.blue[0]-20 <= p['nm'] <= args.blue[1]+20:
            marker = ' 🟡 near blue'
        print(f"   {i+1}. {p['nm']:.1f} nm ({p['ev']:.3f} eV)  I={p['intensity']:.4f}{marker}")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
