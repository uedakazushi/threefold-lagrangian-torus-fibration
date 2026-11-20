#!/usr/bin/env python3
"""Plot the explicit cyclic split-merge multigraph, not an inferred global base."""
from pathlib import Path
import math
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/general_type_matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parent
GENUS = 10
n = GENUS - 1
V = 2 * n
angles = [math.pi / 2 - 2 * math.pi * j / V for j in range(V)]
positions = [(math.cos(t), math.sin(t)) for t in angles]
fig, ax = plt.subplots(figsize=(8.2, 8.3))
fig.patch.set_facecolor('white')
ax.set_aspect('equal')

def edge(i, j, radius, color):
    a = angles[i]
    b = angles[j]
    if j == 0:
        b -= 2 * math.pi
    mid = (a + b) / 2
    start, end = positions[i], positions[j]
    control = (radius * math.cos(mid), radius * math.sin(mid))
    path = MplPath([start, control, end], [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3])
    ax.add_patch(PathPatch(path, facecolor='none', edgecolor=color, lw=2.0))
    t = 0.57
    point = tuple((1-t)**2 * start[k] + 2*t*(1-t)*control[k] + t*t*end[k] for k in range(2))
    deriv = tuple(2*(1-t)*(control[k]-start[k]) + 2*t*(end[k]-control[k]) for k in range(2))
    length = math.hypot(*deriv)
    delta = tuple(0.035*v/length for v in deriv)
    ax.add_patch(FancyArrowPatch(tuple(point[k]-delta[k] for k in range(2)),
                                tuple(point[k]+delta[k] for k in range(2)),
                                arrowstyle='-|>', mutation_scale=11, lw=1.6, color=color))

for k in range(n):
    edge(2*k, 2*k+1, 1.28, '#3e609a')
    edge(2*k, 2*k+1, 0.75, '#3e609a')
    edge(2*k+1, (2*k+2) % V, 1.015, '#5e6972')

for j, (x, y) in enumerate(positions):
    is_split = j % 2 == 0
    ax.scatter([x], [y], s=105, marker='o' if is_split else 's',
               c='#127d83' if is_split else '#c46d31', edgecolors='white', linewidths=1.2, zorder=5)
    label = ('S' if is_split else 'M') + '_{' + str(j//2+1) + '}'
    ax.text(1.23*x, 1.23*y, '$'+label+'$', ha='center', va='center', fontsize=11, color='#25364a')

ax.text(0, 0.18, r'$\Gamma_{10}$', ha='center', va='center', fontsize=29, color='#1b364e')
ax.text(0, -0.04, r'$V=18,\quad E=27$', ha='center', va='center', fontsize=17, color='#1b364e')
ax.text(0, -0.22, r'$b_1=E-V+1=10$', ha='center', va='center', fontsize=16, color='#1b364e')
ax.set_xlim(-1.43, 1.43)
ax.set_ylim(-1.46, 1.43)
ax.axis('off')
ax.set_title('Closed-curve Reeb graph: the sextic case', fontsize=17, pad=18, color='#1b364e')
ax.scatter([-0.70], [-1.38], s=90, c='#127d83', marker='o')
ax.text(-0.62, -1.38, 'Split', va='center', fontsize=12)
ax.scatter([0.06], [-1.38], s=90, c='#c46d31', marker='s')
ax.text(0.14, -1.38, 'Merge', va='center', fontsize=12)
fig.text(0.5, 0.025, '9 split-merge blocks; every vertex is trivalent; no endpoints.',
         ha='center', fontsize=11, color='#46515e')
fig.subplots_adjust(left=0.035, right=0.965, top=0.91, bottom=0.065)
fig.savefig(ROOT/'closed_curve_graph.pdf', bbox_inches='tight')
fig.savefig(ROOT/'closed_curve_graph.png', dpi=180, bbox_inches='tight')
plt.close(fig)
print('Created closed_curve_graph.pdf and closed_curve_graph.png')
