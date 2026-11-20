#!/usr/bin/env python3
"""Reproduce the static research figure; requires matplotlib (checks do not)."""
from pathlib import Path
import math
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'general_type_fibration_mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = Path(__file__).resolve().parent
h, j, c = 4*math.pi, math.pi, 0.3
lo, hi = max(0, j), h/2
u0 = (h+3*j+math.sqrt(h*h-2*h*j+9*j*j))/8
def profile(u):
    return (h-2*u)*math.sqrt(max(0, u*(u-j)))/(math.pi**2)
F = profile(u0)
qs, qm = (F-c)**2, (F+c)**2
assert F > 2*c

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9.0,
    'axes.labelsize': 9.0, 'axes.titlesize': 10.5,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 3.25),
                             gridspec_kw={'width_ratios': [1.25, 1]})
blue, teal, orange = '#235878', '#238A84', '#AD562C'
for a,b,color in [(lo,u0,blue),(u0,hi,teal)]:
    us = [a+(b-a)*i/300 for i in range(301)]
    ax.plot([u/math.pi for u in us], [profile(u) for u in us], color=color, lw=2.3)
ax.axhline(c, color='#697B87', lw=0.9, ls='--')
ax.axvline(u0/math.pi, color='#BAC5CC', lw=0.8, ls=':')
ax.scatter([u0/math.pi], [F], s=25, color=orange, zorder=4)
ax.annotate('$F$', (u0/math.pi,F), xytext=(6,5), textcoords='offset points', color=orange)
ax.text(1.02, c+0.035, '$c$', color='#596B78')
ax.text(1.10, 0.12, 'lower sheet', color=blue)
ax.text(1.62, 0.12, 'upper sheet', color=teal)
ax.set(xlim=(0.98,2.02), ylim=(-0.035,F*1.17),
       xlabel=r'$u/\pi$', ylabel=r'$|W|=f_{h,j}(u)$',
       title='A. Two radial sheets')
ax.set_xticks([1, u0/math.pi, 2], ['$1$', r'$u_0/\pi$', '$2$'])
ax.text(0.5, 1.02, r'$h=4\pi,\ j=\pi,\ c=0.3$',
        transform=ax.transAxes, ha='center', fontsize=8, color='#596B78')
ax.set_title('A. Two radial sheets', pad=23)

bx.plot([-0.60,0], [0,qs], lw=2.3, color=blue)
bx.plot([0.60,0], [0,qs], lw=2.3, color=teal)
bx.plot([0,0], [qs,qm], lw=2.3, color='#425767')
bx.scatter([-0.60,0.60,0], [0,0,qm], s=24, color='#425767', zorder=4)
bx.scatter([0], [qs], s=30, color=orange, zorder=5)
orb_x = 0.60*(1-c*c/qs)
bx.scatter([orb_x], [c*c], s=32, marker='D', color=teal, zorder=5)
bx.annotate('order-2 orbit\nlocal valence 2', (orb_x,c*c),
            xytext=(16,28), textcoords='offset points', fontsize=8,
            arrowprops={'arrowstyle':'-', 'color':teal, 'lw':0.8})
bx.annotate('smooth saddle', (0,qs), xytext=(-84,18),
            textcoords='offset points', fontsize=8, color=orange,
            arrowprops={'arrowstyle':'-', 'color':orange, 'lw':0.8})
bx.text(0.05, qm-0.035, 'maximum', va='top', fontsize=8, color='#425767')
bx.text(-0.60, -0.11, 'minimum', ha='center', fontsize=8, color=blue)
bx.text(0.60, -0.11, 'minimum', ha='center', fontsize=8, color=teal)
bx.set(xlim=(-1.18,1.30), ylim=(-0.16,qm*1.08),
       ylabel=r'$q=Q$', title='B. Connected-level quotient')
bx.set_xticks([])
bx.set_title('B. Connected-level quotient', pad=23)
bx.set_yticks([0, c*c, qs, qm], ['$0$', '$c^2$', '$(F-c)^2$', '$(F+c)^2$'])
bx.tick_params(axis='y', length=0, labelsize=8)
for axis in (ax,bx):
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    for side in ('left','bottom'):
        axis.spines[side].set_color('#9AAAB5')
    axis.tick_params(colors='#425767')
bx.spines['bottom'].set_visible(False)
fig.subplots_adjust(left=0.08, right=0.97, top=0.82, bottom=0.18, wspace=0.49)
fig.savefig(out/'positive_winding_geometry.pdf', bbox_inches='tight', pad_inches=0.06,
            metadata={'Title': 'Positive winding: radial sheets and Reeb branching',
                      'Author': '', 'CreationDate': None, 'ModDate': None})
print('Saved positive_winding_geometry.pdf')
