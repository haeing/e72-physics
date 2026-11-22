"""Crystal Ball 734 MeV/c data from sensitivity/macro/crystalball.cc.

Fit polynomials using the supplied cross-section uncertainties as weights.
Horizontal errors indicate angular bin half-widths and are not fit weights.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    x = np.array([-0.95, -0.80, -0.60, -0.35, 0., 0.35, 0.60, 0.80, 0.95])
    xerr = np.array([0.05, 0.10, 0.10, 0.15, 0.20, 0.15, 0.10, 0.10, 0.05])
    y = np.array([0.17, 0.11, 0.10, 0.04, 0.06, 0.07, 0.09, 0.16, 0.23])
    yerr = np.array([0.05, 0.03, 0.03, 0.01, 0.02, 0.02, 0.03, 0.03, 0.06])

    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color='#333333',
                markersize=5, capsize=3, label='Crystal Ball')
    xx = np.linspace(-1, 1, 500)
    for degree, color, label in [(0, '#335b88', 'Constant fit'),
                                  (2, '#c77932', 'Quadratic fit')]:
        coefficients = np.polyfit(x, y, degree, w=1 / yerr)
        parameter_count = degree + 1
        chi2 = np.sum(((y - np.polyval(coefficients, x)) / yerr) ** 2)
        ndf = len(x) - parameter_count
        print(f'{label}: coefficients (highest power first) = {coefficients}; '
              f'chi2/ndf = {chi2:.3f}/{ndf}')
        ax.plot(xx, np.polyval(coefficients, xx), color=color, lw=2.2, label=label)

    ax.set_xlim(-1, 1)
    ax.set_ylim(0, 0.31)
    ax.set_xticks(np.linspace(-1, 1, 11))
    ax.set_xlabel(r'$\cos\theta_{\eta}^{\mathrm{CM}}$', fontsize=16)
    ax.set_ylabel(r'$d\sigma/d\Omega$ [mb/sr]', fontsize=16)
    ax.tick_params(labelsize=12)
    ax.legend(fontsize=12, frameon=False)
    fig.tight_layout()
    directory = Path(__file__).resolve().parent
    for extension in ('png', 'pdf', 'svg'):
        fig.savefig(directory / f'crystalball_734.{extension}', dpi=240)
    plt.close(fig)


if __name__ == '__main__':
    main()
