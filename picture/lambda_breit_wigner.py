import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Approximate central values from PDG 2024 summary tables.
# Widths are FULL widths Gamma (FWHM), not Gaussian sigma.
resonances = [
    (r'$\Lambda(1520)$', r'$3/2^-$', 1519, 16),
    (r'$\Lambda(1600)$', r'$1/2^+$', 1600, 200),
    (r'$\Lambda(1670)$', r'$1/2^-$', 1674, 30),
    (r'$\Lambda(1690)$', r'$3/2^-$', 1690, 70),
]

# Non-relativistic, constant-width Breit-Wigner (Lorentzian).
# Equal peak amplitudes in separate spin columns; vertical axis is mass.
def breit_wigner(m, m0, gamma):
    return (gamma / 2)**2 / ((m - m0)**2 + (gamma / 2)**2)

fig, ax = plt.subplots(figsize=(9.2, 7.4))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
m = np.linspace(1450, 1800, 2500)
amplitude = 0.42
spin_bases = {'1/2': 0.10, '3/2': 1.25}

for name, jp, m0, gamma in resonances:
    spin = '1/2' if '1/2' in jp else '3/2'
    x_base = spin_bases[spin]
    highlight = '#335b88' if spin == '1/2' else '#c77932'
    intensity = breit_wigner(m, m0, gamma)
    x = x_base + amplitude * intensity
    ax.plot(x, m, lw=2.2, color=highlight)
    ax.fill_betweenx(m, x_base, x, color=highlight, alpha=0.10)
    ax.plot([x_base, x_base + amplitude], [m0, m0], color=highlight,
            alpha=0.30, lw=0.8, ls='--')
    ax.text(x_base + 0.50, m0, name, va='center', fontsize=14,
            color=highlight)

ax.axhline(1663.5, ls=':', lw=1.3, color='#999999', zorder=0)
ax.set_ylim(1480, 1760)
ax.set_xlim(0, 2.50)
ax.set_xticks([base + amplitude / 2 for base in spin_bases.values()],
              [r'$1/2$', r'$3/2$'], fontsize=14)
ax.set_xlabel(r'Spin $J$', fontsize=14)
spin_label_x = np.mean(list(spin_bases.values())) + amplitude / 2
ax.xaxis.set_label_coords(spin_label_x, -0.025, transform=ax.get_xaxis_transform())
ax.set_ylabel(r'Mass (MeV/$c^2$)', fontsize=14)
ax.set_title(r'$\Lambda$ resonances near 1670 MeV', fontsize=17, pad=15)
ax.spines[['top', 'right', 'bottom']].set_visible(False)
ax.spines['left'].set_position(('data', 0.10))
ax.tick_params(axis='x', length=0, pad=10)
ax.tick_params(axis='y', labelsize=11)
fig.subplots_adjust(left=0.15, right=0.97, top=0.91, bottom=0.12)
output_dir = Path(__file__).resolve().parent
fig.savefig(output_dir / 'lambda_breit_wigner.png', dpi=240, bbox_inches='tight')
fig.savefig(output_dir / 'lambda_breit_wigner.svg', bbox_inches='tight')
plt.show()
