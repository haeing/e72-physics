"""Lambda/Sigma resonance level chart using PDG 2025 estimates.

Run: python pdg_hyperon_resonances.py [--no-show]
Vertical bars span M - Gamma/2 to M + Gamma/2; the horizontal mark is M.
Select states by representative mass (1500--1800 MeV), one bar per multiplet.
Arrows indicate widths continuing beyond the visible energy range.
PDG estimate ranges in CSV are not statistical errors. Width is the full
resonance width, not mass uncertainty. Sources are an offline 2025 snapshot.
Writes PNG, SVG and CSV. Requires numpy and matplotlib.
"""
import argparse
import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

BASE = 'https://pdg.lbl.gov/2025/'
LAMBDA_SOURCE = BASE + 'tables/rpp2025-tab-baryons-Lambda.pdf'
SIGMA_SOURCE = BASE + 'tables/rpp2025-tab-baryons-Sigma.pdf'
SIGMA1620_SOURCE = BASE + 'listings/rpp2025-list-sigma-1620.pdf'


@dataclass(frozen=True)
class Resonance:
    family: str
    nominal: int
    jp: str
    mass: float
    width: float
    mass_low: float
    mass_high: float
    width_low: float
    width_high: float
    source: str
    note: str = "PDG estimate ranges"


# All dimensional values are MeV. Use PDG's approximate value, not necessarily
# the midpoint of the quoted range or the nominal mass in the particle name.
RESONANCES = (
    Resonance('Lambda', 1520, '3/2-', 1519, 16, 1518, 1520, 15, 17, LAMBDA_SOURCE),
    Resonance('Lambda', 1600, '1/2+', 1600, 200, 1570, 1630, 150, 250, LAMBDA_SOURCE),
    Resonance('Lambda', 1670, '1/2-', 1674, 30, 1670, 1678, 25, 35, LAMBDA_SOURCE),
    Resonance('Lambda', 1690, '3/2-', 1690, 70, 1685, 1695, 60, 80, LAMBDA_SOURCE),
    Resonance('Lambda', 1710, '1/2+', 1713, 180, 1700, 1726, 138, 222,
              BASE + 'listings/rpp2025-list-lambda-1710.pdf', 'ZHANG 13A single measurement; bounds are quoted +/- errors'),
    Resonance('Lambda', 1810, '1/2+', 1790, 110, 1740, 1840, 50, 170, LAMBDA_SOURCE),
    Resonance('Lambda', 1800, '1/2-', 1800, 200, 1750, 1850, 150, 250, LAMBDA_SOURCE),
    Resonance('Sigma', 1580, '3/2-', 1580, 11, 1580, 1580, 7, 15,
              BASE + 'listings/rpp2025-list-sigma-1580.pdf', 'Mass: PDG estimate; width: LITCHFIELD 74, 11 +/- 4 MeV; no PDG recommended width'),
    Resonance('Sigma', 1620, '1/2-', 1620, 70, 1600, 1650, 40, 100, SIGMA1620_SOURCE),
    Resonance('Sigma', 1660, '1/2+', 1660, 200, 1640, 1680, 100, 300, SIGMA_SOURCE),
    Resonance('Sigma', 1670, '3/2-', 1675, 70, 1665, 1685, 40, 100, SIGMA_SOURCE),
    Resonance('Sigma', 1750, '1/2-', 1750, 150, 1700, 1800, 100, 200, SIGMA_SOURCE),
    Resonance('Sigma', 1775, '5/2-', 1775, 120, 1770, 1780, 105, 135, SIGMA_SOURCE),
    Resonance('Sigma', 1780, '3/2+', 1780, 200, 1730, 1830, 100, 300,
              BASE + 'listings/rpp2025-list-sigma-1780-3o2-plus.pdf'),
)


def breit_wigner(energy, mass, width):
    """Peak = 1, full width at half maximum = width."""
    half = width / 2
    return half**2 / ((np.asarray(energy) - mass)**2 + half**2)


def bin_averages(edges, mass, width):
    half = width / 2
    primitive = half * np.arctan((edges - mass) / half)
    return np.diff(primitive) / np.diff(edges)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ymin', type=float, default=1500, help='Minimum CM energy [MeV]')
    parser.add_argument('--ymax', type=float, default=1800, help='Maximum CM energy [MeV]')
    parser.add_argument('--output', default='pdg_hyperon_resonances', help='Output filename stem')
    parser.add_argument('--no-show', action='store_true', help='Save without opening a window')
    args = parser.parse_args()
    if not all(np.isfinite(v) for v in (args.ymin, args.ymax)) or args.ymax <= args.ymin:
        parser.error('Require finite energy limits with ymax > ymin')

    import matplotlib
    if args.no_show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    selected = [r for r in RESONANCES if args.ymin <= r.mass <= args.ymax]
    if not selected:
        parser.error('No stored resonances in the requested energy interval')
    fig, ax = plt.subplots(figsize=(8, 12))
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.15, top=0.92)
    bar_width = 0.40
    family_centers = []
    x = 0
    divider = None
    for family, cmap in [('Lambda', plt.cm.Blues), ('Sigma', plt.cm.Oranges)]:
        members = [r for r in selected if r.family == family]
        start = x
        for i, r in enumerate(members):
            color = cmap(0.5 + 0.4*i/max(1, len(members)-1))
            low, high = r.mass-r.width/2, r.mass+r.width/2
            ax.bar(x, r.width, bottom=low, width=bar_width,
                   facecolor=color, edgecolor=color, alpha=0.4, linewidth=1, zorder=2)
            for edge in (low, high):
                if args.ymin <= edge <= args.ymax:
                    ax.hlines(edge, x-bar_width/2, x+bar_width/2, color=color, linewidth=1.5, zorder=3)
            ax.hlines(r.mass, x-bar_width/2, x+bar_width/2, color=color, linewidth=2, zorder=3)
            for bound, end, marker in [(args.ymin, low, 'v'), (args.ymax, high, '^')]:
                if (marker == 'v' and end < bound) or (marker == '^' and end > bound):
                    ax.plot(x, bound, marker=marker, color=color, markersize=5, clip_on=False, zorder=4)
            symbol = r'\Lambda' if family == 'Lambda' else r'\Sigma'
            spin, parity = r.jp[:-1], r.jp[-1]
            dagger = '†' if (family, r.nominal) in [('Lambda', 1710), ('Sigma', 1580)] else ''
            label = (rf'${symbol}({r.nominal})\;{spin}^{{{parity}}}$' + dagger
                     + f'\nM={r.mass:g}, Γ={r.width:g}')
            # Put labels beside the bars, within their own narrow lanes.
            text_y = np.clip(r.mass, args.ymin+37, args.ymax-37)
            ax.text(x+0.29, text_y, label, rotation=90, ha='center', va='center',
                    fontsize=8, color='0.15', zorder=5)
            x += 1
        if members:
            family_centers.append(((start+x-1)/2, r'$\Lambda$' if family == 'Lambda' else r'$\Sigma$'))
        if family == 'Lambda':
            divider = x-0.1
            x += 0.8
    if divider is not None:
        ax.axvline(divider, color='0.8', linewidth=1)
    ax.set(xlim=(-0.55, x-0.3), ylim=(args.ymin, args.ymax), ylabel='CM energy W [MeV]')
    ax.set_title(r'$\Lambda$ and $\Sigma$ resonance masses and widths', fontsize=14, pad=18)
    ax.set_xticks([p[0] for p in family_centers], [p[1] for p in family_centers], fontsize=22)
    ax.tick_params(axis='x', length=0, pad=12)
    ax.set_yticks(np.arange(np.ceil(args.ymin/50)*50, args.ymax+1, 50))
    ax.grid(axis='y', alpha=0.2, zorder=0)
    ax.spines[['top', 'right']].set_visible(False)
    fig.text(0.12, 0.085, r'Bar: $M\pm\Gamma/2$; line: $M$; triangles: width continues off-scale.', fontsize=9)
    fig.text(0.12, 0.059, 'PDG 2025 • representative masses 1500–1800 MeV; widths in MeV.', fontsize=9)
    fig.text(0.12, 0.033, '† Individual measurements: Λ(1710), ZHANG 13A;\n   Σ(1580) width, LITCHFIELD 74 (no PDG recommended width).', fontsize=8, color='0.35')
    stem = Path(args.output)
    for suffix in ('.png', '.svg'):
        output = stem.with_suffix(suffix)
        fig.savefig(output, dpi=180)
        print(f'Saved: {output}')
    output = stem.with_suffix('.csv')
    with output.open('w', newline='') as file:
        fields = list(asdict(RESONANCES[0]))
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(asdict(r) for r in selected)
    print(f'Saved: {output} (all mass/width columns in MeV)')
    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
