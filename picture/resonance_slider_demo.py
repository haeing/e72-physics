"""Interactive two-resonance partial-wave toy model.

Run: python resonance_slider_demo.py
Save automatic slider sweep: python resonance_slider_demo.py --save resonance_demo.mp4
MP4 export requires ffmpeg; --seconds and --fps control duration and frame rate.
Requires numpy, matplotlib and an interactive GUI backend (TkAgg/QtAgg).
Illustrative constant-width amplitudes, not a fit to K-p -> Lambda eta data.
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons

E = np.linspace(1.58, 1.82, 260)
cos_theta = np.linspace(-1, 1, 250)
E0 = 1.700
# Spin-0 + spin-1/2 scattering: labels map to (orbital l, J=l +/- 1/2).
# Parity is (-1)^(l+1) for a pseudoscalar meson and positive-parity baryon.
WAVES = {'1/2-': (0, 'plus'), '1/2+': (1, 'minus'),
         '3/2-': (2, 'minus'), '3/2+': (1, 'plus')}
COLORS = ('tab:blue', 'tab:orange')


def bw(energy, mass, width):
    return (width / 2) / (mass - energy - 0.5j * width)


fig = plt.figure(figsize=(15, 9), layout=None)
fig.suptitle('Two-resonance toy model: cross sections and recoil polarization', fontsize=16, y=0.975)
ax_tot = fig.add_axes([0.40, 0.59, 0.235, 0.30])
ax_diff = fig.add_axes([0.735, 0.59, 0.235, 0.30])
ax_pol = fig.add_axes([0.735, 0.14, 0.235, 0.30])
slider_specs = [
    ('M1 [GeV]', 1.60, 1.78, 1.670),
    ('Γ1 [GeV]', 0.015, 0.16, 0.055),
    ('M2 [GeV]', 1.60, 1.78, 1.720),
    ('Γ2 [GeV]', 0.015, 0.16, 0.095),
    ('Strength 2', 0.0, 2.0, 0.8),
    ('Phase 2 [deg]', -180, 180, 75),
    ('Energy [GeV]', 1.59, 1.81, E0),
]
sliders = {}
for i, (label, lo, hi, initial) in enumerate(slider_specs):
    y = 0.85 - 0.08*i
    ax = fig.add_axes([0.055, y, 0.23, 0.022])
    slider = Slider(ax, '', lo, hi, valinit=initial,
                    valfmt='%1.0f' if 'deg' in label else '%1.3f')
    fig.text(0.055, y + 0.032, label, fontsize=10)
    sliders[label] = slider

radio1_ax = fig.add_axes([0.045, 0.08, 0.125, 0.19])
radio2_ax = fig.add_axes([0.195, 0.08, 0.125, 0.19])
radio1_ax.set_title(r'Res. 1 $J^P$', fontsize=11, color=COLORS[0])
radio2_ax.set_title(r'Res. 2 $J^P$', fontsize=11, color=COLORS[1])
radio1 = RadioButtons(radio1_ax, tuple(WAVES), active=0)
radio2 = RadioButtons(radio2_ax, tuple(WAVES), active=1)


def plot_lines(ax):
    total, = ax.plot([], [], lw=2.5, color='black', label='Combined')
    parts = [ax.plot([], [], lw=1.8, ls='--', color=color,
                     label=f'Res. {i + 1} only')[0] for i, color in enumerate(COLORS)]
    ax.legend(loc='upper right', fontsize=8)
    return total, parts


line_total, total_parts = plot_lines(ax_tot)
line_energy = ax_tot.axvline(E0, color='gray', ls=':', alpha=0.7)
line_diff, diff_parts = plot_lines(ax_diff)
line_pol, pol_parts = plot_lines(ax_pol)
ax_tot.set(xlabel='CM energy W [GeV]', ylabel='Total cross section [arb. units]', title='Total cross section')
ax_diff.set(xlabel=r'$\cos\theta$', ylabel=r'$d\sigma/d\Omega$ [arb. units]', title='Differential cross section')
ax_pol.set(xlabel=r'$\cos\theta$', ylabel=r'$P_\Lambda$', title='Recoil polarization', ylim=(-1.05, 1.05))
ax_pol.axhline(0, color='gray', lw=0.8)
ax_diff.set_xlim(-1, 1)
ax_pol.set_xlim(-1, 1)


def resonance_amplitudes(energy):
    a1 = bw(energy, sliders['M1 [GeV]'].val, sliders['Γ1 [GeV]'].val)
    a2 = (sliders['Strength 2'].val * np.exp(1j * np.deg2rad(sliders['Phase 2 [deg]'].val))
          * bw(energy, sliders['M2 [GeV]'].val, sliders['Γ2 [GeV]'].val))
    return a1, a2


def angular_amplitudes(amplitude, label, x):
    """f=sum[(l+1)t_l+ + l*t_l-]P_l; g=sum[t_l- - t_l+]sin(theta)P'_l."""
    l, branch = WAVES[label]
    polynomial = (np.ones_like(x), x, (3*x*x - 1)/2)[l]
    derivative = (np.zeros_like(x), np.ones_like(x), 3*x)[l]
    plus = branch == 'plus'
    f = (l + 1 if plus else l) * amplitude * polynomial
    g = (-1 if plus else 1) * amplitude * np.sqrt(np.maximum(0, 1-x*x)) * derivative
    return f, g


def integrated_cross_section(amplitudes, labels):
    # Equal partial waves interfere; different (l,J) waves integrate orthogonally.
    grouped = {}
    for amplitude, label in zip(amplitudes, labels):
        grouped[label] = grouped.get(label, 0) + amplitude
    return 4*np.pi*sum((l + 1 if branch == 'plus' else l)*np.abs(grouped[label])**2
                       for label, (l, branch) in WAVES.items() if label in grouped)


def observables(f, g):
    intensity = np.abs(f)**2 + np.abs(g)**2
    polarization = np.divide(2*np.imag(f*np.conj(g)), intensity,
                             out=np.zeros_like(intensity), where=intensity > 1e-12)
    return intensity, polarization


def update(_=None):
    labels = (radio1.value_selected, radio2.value_selected)
    amplitudes = resonance_amplitudes(E)
    totals = [integrated_cross_section(amplitudes, labels)]
    totals.extend(integrated_cross_section([a], [label]) for a, label in zip(amplitudes, labels))
    for line, values in zip([line_total, *total_parts], totals):
        line.set_data(E, values)
    ax_tot.set_xlim(E.min(), E.max())
    ax_tot.set_ylim(0, max(0.1, 1.15*max(np.max(v) for v in totals)))
    energy = sliders['Energy [GeV]'].val
    line_energy.set_xdata([energy, energy])
    selected = resonance_amplitudes(energy)
    components = [angular_amplitudes(a, label, cos_theta) for a, label in zip(selected, labels)]
    combined = tuple(sum(pair[i] for pair in components) for i in range(2))
    results = [observables(*pair) for pair in [combined, *components]]
    for dline, pline, (intensity, polarization) in zip(
            [line_diff, *diff_parts], [line_pol, *pol_parts], results):
        dline.set_data(cos_theta, intensity)
        pline.set_data(cos_theta, polarization)
    ax_diff.set_ylim(0, max(0.1, 1.15*max(np.max(r[0]) for r in results)))
    fig.canvas.draw_idle()


for slider in sliders.values():
    slider.on_changed(update)
radio1.on_clicked(update)
radio2.on_clicked(update)
update()
def save_animation(filename, seconds=14, fps=15):
    """Sweep each slider in turn: initial -> minimum -> maximum -> initial.

    Widget positions and values are rendered along with the three plots.
    Spin selections stay fixed; all sliders are restored after export.
    """
    from matplotlib.animation import FFMpegWriter

    if seconds <= 0 or fps <= 0:
        raise ValueError('seconds and fps must be positive')
    if not FFMpegWriter.isAvailable():
        raise RuntimeError('MP4 export requires ffmpeg on PATH.')
    if Path(filename).suffix.lower() != '.mp4':
        raise ValueError('Use an .mp4 output filename.')
    controls = list(sliders.values())
    originals = [slider.val for slider in controls]
    draw_states = [slider.drawon for slider in controls]
    event_states = [slider.eventson for slider in controls]
    frames = max(2, round(seconds * fps))
    writer = FFMpegWriter(fps=fps, codec='libx264',
                         extra_args=['-pix_fmt', 'yuv420p', '-crf', '20'])
    try:
        for slider in controls:
            slider.drawon = False
            slider.eventson = False
        with writer.saving(fig, str(filename), dpi=100):
            for frame in range(frames):
                progress = frame / (frames - 1) * len(controls)
                index = min(int(progress), len(controls) - 1)
                local = min(progress - index, 1.0)
                slider = controls[index]
                stops = (originals[index], slider.valmin, slider.valmax, originals[index])
                stage = min(int(local * 3), 2)
                fraction = (local * 3 - stage)
                blend = (1 - np.cos(np.pi * fraction)) / 2
                for control, original in zip(controls, originals):
                    control.set_val(original)
                slider.set_val(stops[stage] + blend * (stops[stage + 1] - stops[stage]))
                update()
                writer.grab_frame()
                if frame % max(1, fps * 2) == 0:
                    print(f'Export: {frame + 1}/{frames} frames', flush=True)
    finally:
        for slider, value, drawon, eventson in zip(controls, originals, draw_states, event_states):
            slider.set_val(value)
            slider.drawon = drawon
            slider.eventson = eventson
        update()
    print(f'Saved: {filename}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--save', metavar='FILE.mp4', help='Export automatic slider sweeps without opening a window')
    parser.add_argument('--seconds', type=float, default=14, help='Video duration (default: 14)')
    parser.add_argument('--fps', type=int, default=15, help='Video frame rate (default: 15)')
    parser.add_argument('--spin1', choices=tuple(WAVES), default='1/2-')
    parser.add_argument('--spin2', choices=tuple(WAVES), default='1/2+')
    args = parser.parse_args()
    if args.seconds <= 0 or args.fps <= 0:
        parser.error('--seconds and --fps must be positive')
    radio1.set_active(tuple(WAVES).index(args.spin1))
    radio2.set_active(tuple(WAVES).index(args.spin2))
    if args.save:
        # A noninteractive canvas lets exports run without a desktop window.
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        FigureCanvasAgg(fig)
        save_animation(args.save, args.seconds, args.fps)
    else:
        plt.show()


if __name__ == '__main__':
    main()
