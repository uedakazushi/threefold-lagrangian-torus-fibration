#!/usr/bin/env python3
"""Scientific illustration of whole-fiber chambers and annular truncation."""
from pathlib import Path
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def root(h, radius, upper=False):
    a, b = ((2*h/3, h) if upper else (0, 2*h/3))
    for _ in range(100):
        m = (a+b)/2
        f = m * math.sqrt(max(0, h-m)) / math.pi**1.5
        if (f < radius) != upper:
            a = m
        else:
            b = m
    return (a+b)/2


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    h, e1, c = 3.0, 1.0, 0.02
    wall = (h-e1)/2
    F = 2*h**1.5/(3*math.sqrt(3)*math.pi**1.5)
    qs = (F-c)**2
    q0 = (c*c+qs)/2
    qmin = c*c/2
    qmax = qs+0.006
    focus_max = root(h, 2*c)
    saddle_min = root(h, F-2*c)
    saddle_max = root(h, F-2*c, True)
    fig, (ax, graph) = plt.subplots(1, 2, figsize=(11.5, 4.6),
                                  gridspec_kw={"width_ratios": [1.2, 1]})
    u = np.linspace(0, h, 900)
    f = u*np.sqrt(h-u)/math.pi**1.5
    ax.axvspan(0, wall, color="#e7f3ee", zorder=0)
    ax.axvspan(wall, h, color="#edf2fa", zorder=0)
    ax.axvspan(0, focus_max, color="#278267", alpha=.4, zorder=1)
    ax.axvspan(saddle_min, saddle_max, color="#df9637", alpha=.33, zorder=1)
    ax.plot(u, f, color="#18374a", lw=2.2)
    ax.axvline(wall, color="#657482", ls="--", lw=1.3)
    ax.axhline(2*c, color="#278267", ls=":", lw=1)
    ax.axhline(F-2*c, color="#b56a13", ls=":", lw=1)
    ax.text(.48, .35, r"$\mathcal{C}_j$", ha="center", color="#276c57")
    ax.text(2.05, .39, r"$\mathcal{C}_i$", ha="center", color="#345e91")
    ax.text(wall+.03, .1, r"$u_b=(h-e_1)/2$", rotation=90, va="bottom", fontsize=9)
    ax.annotate("whole focus", xy=(focus_max/2, .022), xytext=(.17, .16),
                arrowprops={"arrowstyle": "-", "color": "#278267"}, color="#216b55")
    ax.annotate("whole saddle", xy=(2, .31), xytext=(2.12, .20),
                arrowprops={"arrowstyle": "-", "color": "#b56a13"}, color="#9c5911")
    ax.set(xlim=(0, h), ylim=(0, .42), xlabel=r"$u=P_y$",
           ylabel=r"$f_h(u)=|W|$", title="A. Whole fibers occupy different chambers")
    ax.text(.02, .97, r"$h=3,\ e_1=1,\ c=0.02$", transform=ax.transAxes,
            va="top", fontsize=9)
    # A fixed (h,j=0) slice of the quotient. Only the lower component
    # is retained below q0; both components are retained above q0.
    # Heights encode order only, to separate the focus from the nearby lower cut.
    height = lambda q: float(np.interp(q, [qmin, c*c, q0, qs, qmax], [0, .2, .48, .85, 1]))
    xleft = lambda q: -(1-height(q)/.85)
    xright = lambda q: 1-height(q)/.85
    graph.plot([xleft(qmin), 0], [height(qmin), height(qs)], color="#276b79", lw=3)
    graph.plot([xright(q0), 0], [height(q0), height(qs)], color="#276b79", lw=3)
    graph.plot([0, 0], [height(qs), height(qmax)], color="#276b79", lw=3)
    graph.plot([xright(qmin), xright(q0)], [height(qmin), height(q0)],
               color="#b5bcc2", lw=1.7, ls="--")
    graph.scatter([xleft(c*c)], [height(c*c)], marker="*", s=120, color="#ae3f44", zorder=5)
    graph.scatter([xright(c*c)], [height(c*c)], marker="o", s=35,
                  facecolors="white", edgecolors="#929ba4", zorder=5)
    graph.scatter([0], [height(qs)], s=55, color="#c1842e", zorder=5)
    graph.scatter([xleft(qmin), xright(q0), 0], [height(qmin), height(q0), height(qmax)],
                  s=28, color="#276b79", zorder=5)
    graph.annotate("compact focus component", (xleft(c*c), height(c*c)),
                   xytext=(-1.18, .34), fontsize=9, color="#ae3f44",
                   arrowprops={"arrowstyle": "-", "color": "#ae3f44"})
    graph.annotate(r"excluded $x=0$ component", (xright(c*c), height(c*c)),
                   xytext=(.23, .08), fontsize=9, color="#777f87",
                   arrowprops={"arrowstyle": "-", "color": "#929ba4"})
    graph.annotate(r"regular cut at $q_0$", (xright(q0), height(q0)),
                   xytext=(.36, .59), fontsize=9,
                   arrowprops={"arrowstyle": "-", "color": "#276b79"})
    graph.text(.1, height(qs)-.02, r"$T^2\times(S^1\vee S^1)$", va="top", fontsize=9)
    graph.set(xlim=(-1.35, 1.35), ylim=(-.05, 1.08),
              yticks=[0, .2, .48, .85, 1],
              yticklabels=[r"$q_-$", r"$c^2$", r"$q_0$", r"$q_s$", r"$q_+$"],
              xticks=[], ylabel="Level order", title="B. Connected-fiber quotient (schematic)")
    graph.spines["bottom"].set_visible(False)
    fig.tight_layout(w_pad=3)
    path=Path(__file__).with_name("positive_transport_geometry.pdf")
    fig.savefig(path, bbox_inches="tight", metadata={"Title": "Positive-winding transport and annular insertion",
                                                     "Author": "", "CreationDate": None, "ModDate": None})
    plt.close(fig)
    print(path)


if __name__ == "__main__":
    main()
