"""Plot Crystal Ball polarization and D03/P03 coordinates from polar_th.cc.

Theory coordinates are connected by natural cubic splines (not fitted to data).
Experimental coordinates, including the final x=0.85, follow the source verbatim.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def smooth(x, y, xx):
    x, y = np.asarray(x), np.asarray(y)
    h = np.diff(x)
    matrix = np.eye(len(x))
    rhs = np.zeros(len(x))
    for i in range(1, len(x) - 1):
        matrix[i, i - 1:i + 2] = [h[i - 1], 2 * (h[i - 1] + h[i]), h[i]]
        rhs[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
    second = np.linalg.solve(matrix, rhs)
    i = np.clip(np.searchsorted(x, xx) - 1, 0, len(x) - 2)
    a, b = (x[i + 1] - xx) / h[i], (xx - x[i]) / h[i]
    return (a * y[i] + b * y[i + 1]
            + ((a**3 - a) * second[i] + (b**3 - b) * second[i + 1]) * h[i]**2 / 6)


def main():
    x = [-0.95, -0.80, -0.60, -0.35, 0., 0.35, 0.60, 0.80, 0.85]
    xerr = [0.05, 0.10, 0.10, 0.15, 0.20, 0.15, 0.10, 0.10, 0.05]
    datasets = {
        735: (
            [0.16, 0.27, -0.26, -0.17, 0.05, 0.05, 0.21, 0.16, -0.02],
            [0.29, 0.23, 0.27, 0.22, 0.22, 0.24, 0.23, 0.23, 0.27],
            ([-1., -0.4797, -0.21805, 0.1594, 0.48271, 0.73383, 1.],
             [0., -0.22509, -0.14823, 0.20741, 0.3371, 0.27142, 0.]),
            ([-1., -0.55779, 0.10178, 1.], [0., -0.64754, -0.56467, 0.])),
        765: (
            [0.06, -0.16, 0.14, -0.02, -0.04, 0.59, 0.07, -0.48, -0.13],
            [0.32, 0.23, 0.26, 0.23, 0.24, 0.26, 0.26, 0.26, 0.38],
            ([-1., -0.80968, -0.48915, -0.05509, 0.28548, 0.62604, 0.74624, 0.84975, 0.95659, 1.],
             [0., -0.03167, -0.0362, 0.02715, 0.09955, 0.14932, 0.1448, 0.1267, 0.0724, 0.]),
            ([-1., -0.43239, -0.00167, 0.28214, 0.69616, 1.],
             [0., -0.38017, -0.53157, -0.47747, -0.23689, 0.]))
    }
    directory = Path(__file__).resolve().parent
    xx = np.linspace(-1, 1, 500)
    for momentum, (y, yerr, d03, p03) in datasets.items():
        fig, ax = plt.subplots(figsize=(8, 6.5))
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color='#333333',
                    markersize=5, capsize=3, label='Crystal Ball')
        ax.axhline(0, color='#377eb8', lw=2.2, linestyle='--',
                   label=r'$\Lambda(1670)$ only', zorder=1)
        for (tx, ty), color, label in [
                (d03, '#9467bd', r'$\Lambda(1670)$ + $J^P = 3/2^-$ resonance'),
                (p03, '#2a9d70', r'$\Lambda(1670)$ + $J^P = 3/2^+$ resonance')]:
            ax.plot(xx, smooth(tx, ty, xx), color=color, lw=2.2, label=label)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_xticks(np.linspace(-1, 1, 11))
        ax.set_xlabel(r'$\cos\theta_{\eta}^{\mathrm{CM}}$', fontsize=16)
        ax.set_ylabel(r'$P_{\Lambda}$', fontsize=16)
        ax.tick_params(labelsize=12)
        ax.legend(fontsize=12, frameon=False, loc='upper left')
        fig.tight_layout()
        for extension in ('png', 'pdf', 'svg'):
            fig.savefig(directory / f'polarization_{momentum}.{extension}', dpi=240)
        plt.close(fig)
        print(f'Saved polarization_{momentum}: 9 experimental points, '
              'Lambda(1670) baseline and additional D03/P03 resonance curves')


if __name__ == '__main__':
    main()
