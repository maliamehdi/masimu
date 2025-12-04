#!/usr/bin/env python3
"""Plot efficiency vs neutron energy (semi-log X) from Tetra summary .txt files.

Usage:
  python3 plot_efficiency.py --dir mono --out mono/mono_eff_scan_py

The script looks for files matching: output_mono_n_*_Tetra_summary.txt
It extracts energy from the filename ("100keV" -> 0.1 MeV, "11" -> 11 MeV),
parses the line "Mean measured multiplicity: <val>", computes efficiency = val*100,
and falls back to the "Events (all detections <= gate): N (eff = X %)" line or
to events/total if needed.
"""

import re
import glob
import os
import argparse
import csv
import math
# import matplotlib only if plotting is requested (avoids hard dependency)

def parse_energy_from_name(fname):
    m = re.search(r"output_mono_n_([^_]+)_Tetra_summary\.txt", fname)
    if not m:
        return None
    token = m.group(1)
    if token.endswith('keV'):
        try:
            kev = float(token[:-3])
            return kev/1000.0
        except Exception:
            return None
    # pure integer -> MeV
    if re.fullmatch(r"\d+", token):
        return float(token)
    try:
        return float(token)
    except Exception:
        return None

def parse_summary(path):
    out = {'file': path, 'E_MeV': None, 'eff_percent': None, 'total': None, 'eventsAllDet': None}
    fname = os.path.basename(path)
    out['E_MeV'] = parse_energy_from_name(fname)
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            m = re.search(r"Mean measured multiplicity:\s*([0-9eE+\.-]+)", line)
            if m:
                try:
                    out['eff_percent'] = float(m.group(1)) * 100.0
                except Exception:
                    pass
            m = re.search(r"^Total events:\s*(\d+)", line)
            if m:
                out['total'] = int(m.group(1))
            m = re.search(r"^Events \(all detections .*\):\s*(\d+) \(eff = ([0-9eE+\.-]+) %\)", line)
            if m:
                out['eventsAllDet'] = int(m.group(1))
                try:
                    if out['eff_percent'] is None:
                        out['eff_percent'] = float(m.group(2))
                except Exception:
                    pass
    # fallback compute
    if out['eff_percent'] is None and out['eventsAllDet'] is not None and out['total']:
        out['eff_percent'] = 100.0 * out['eventsAllDet'] / out['total']
    return out

def main(directory, outprefix, do_plot=True):
    pattern = os.path.join(directory, 'output_mono_n_*_Tetra_summary.txt')
    files = sorted(glob.glob(pattern))
    if not files:
        print('No summary files found with pattern', pattern)
        return 1
    points = []
    for p in files:
        info = parse_summary(p)
        if info['E_MeV'] is None:
            print('Skipping (no energy parsed):', p)
            continue
        if info['eff_percent'] is None:
            print('Skipping (no efficiency parsed):', p)
            continue
        # compute uncertainty: sqrt(meanMeasuredMultiplicity / totalEvents) * 100
        eff_unc = None
        if info.get('total') and info.get('total')>0 and 'Mean measured multiplicity' in open(p).read():
            # try to extract mean measured multiplicity directly
            mval = None
            with open(p,'r') as fh:
                for L in fh:
                    mm = re.search(r"Mean measured multiplicity:\s*([0-9eE+\.-]+)", L)
                    if mm:
                        mval = float(mm.group(1)); break
            if mval is not None:
                eff_unc = math.sqrt(mval / float(info['total'])) * 100.0
        # fallback: use binomial uncertainty on fraction eventsAllDet/total
        if eff_unc is None and info.get('eventsAllDet') and info.get('total') and info['total']>0:
            frac = float(info['eventsAllDet'])/float(info['total'])
            eff_unc = math.sqrt(frac*(1.0-frac)/float(info['total'])) * 100.0
        info['eff_unc_percent'] = eff_unc if eff_unc is not None else 0.0
        points.append(info)
    if not points:
        print('No valid points to plot')
        return 1
    points.sort(key=lambda x: x['E_MeV'])

    energies = [p['E_MeV'] for p in points]
    effs = [p['eff_percent'] for p in points]
    errs = [p.get('eff_unc_percent',0.0) for p in points]

    if do_plot:
        try:
            import matplotlib.pyplot as plt
        except Exception as e:
            print('matplotlib not available, skipping plot:', e)
            do_plot = False

    if do_plot:
        # Plot semi-log X with error bars
        plt.figure(figsize=(10,6))
        plt.errorbar(energies, effs, yerr=errs, marker='o', linestyle='-', capsize=3)
        plt.xscale('log')
        plt.grid(which='both', linestyle='--', linewidth=0.5)
        plt.xlabel('Neutron energy (MeV)')
        plt.ylabel('Detection efficiency (%)')
        plt.title('Detection efficiency (Mean measured multiplicity *100)')
        # Save
        png = outprefix + '_eff_vs_energy_py.png'
        plt.savefig(png, dpi=200)
    # also write CSV
    csvf = outprefix + '_eff_vs_energy_py.csv'
    with open(csvf, 'w', newline='') as cf:
        writer = csv.writer(cf)
        writer.writerow(['E_MeV', 'eff_percent', 'eff_unc_percent', 'totalEvents', 'eventsAllDet', 'file'])
        for p in points:
            writer.writerow([p['E_MeV'], p['eff_percent'], p.get('eff_unc_percent',0.0), p.get('total',''), p.get('eventsAllDet',''), p['file']])
    if do_plot:
        print('Wrote', png, 'and', csvf)
    else:
        print('Wrote CSV:', csvf, '(plot skipped)')
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot efficiency vs neutron energy (semi-log X)')
    parser.add_argument('--dir', default='mono', help='Directory containing summary .txt files')
    parser.add_argument('--out', default='mono/mono_eff_scan_py', help='Output prefix (no extension)')
    parser.add_argument('--no-plot', action='store_true', help='Do not attempt to import matplotlib or draw plot (only write CSV)')
    args = parser.parse_args()
    do_plot = not args.no_plot
    raise SystemExit(main(args.dir, args.out, do_plot=do_plot))
