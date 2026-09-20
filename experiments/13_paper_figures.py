"""Figures for the paper, generated from the run outputs in results/.

fig1  the confound and the fix (data: results/figures.json, exps 01-02)
fig2  design guidance: power vs plots; stacked vs averaged repeat surveys
      (data: results/06_power_pegasus.txt, results/11_design_extras.txt,
      transcribed below with their settings)
fig3  the Everglades application: standardized chi window curves along the
      TP axis, with z-hat, TITAN's change points, and the 10 ug/L criterion
      (window curves recomputed here exactly as in 09_glades.py)

Colors are the project identity (copper / patina / warm gray); every
multi-series panel also separates series by line style and direct labels,
so identity survives grayscale print and color-vision deficiency.
Output: paper/figs/*.pdf (+ .png previews).
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COPPER, PATINA, GRAY, INK = "#A6552C", "#3E7C74", "#6E675F", "#23201C"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGS = os.path.join(ROOT, "paper", "figs")
os.makedirs(FIGS, exist_ok=True)

plt.rcParams.update({
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.edgecolor": GRAY, "axes.labelcolor": INK,
    "xtick.color": GRAY, "ytick.color": GRAY,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#DDD8D2", "grid.linewidth": 0.5,
    "legend.frameon": False,
})


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print(name)


# ---- fig 1: the confound and the fix ---------------------------------------
fj = json.load(open(os.path.join(ROOT, "results", "figures.json")))
rs, md = fj["radius_sweep"], fj["matched_density"]

fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.5))
a.plot(rs["r"], rs["d_density"], color=COPPER, lw=1.8, ls="-")
a.plot(rs["r"], rs["d_chi"], color=PATINA, lw=1.8, ls="--")
a.plot(rs["r"], rs["d_dev"], color=GRAY, lw=1.6, ls="-.")
a.annotate("edge density", (rs["r"][3], rs["d_density"][3]),
           xytext=(0, 6), textcoords="offset points", color=INK, fontsize=8)
a.annotate(r"$\chi$ (flag complex)", (rs["r"][10], rs["d_chi"][10]),
           xytext=(0, 6), textcoords="offset points", color=INK, fontsize=8)
a.annotate(r"$\chi-$ER expectation", (rs["r"][6], rs["d_dev"][6]),
           xytext=(2, -13), textcoords="offset points", color=INK, fontsize=8)
a.set_xlabel("similarity radius $r$")
a.set_ylabel(r"$|$Cohen $d|$ between groups")
a.set_title("(a)  similarity-indexed: density dominates", loc="left")

b.plot(md["m"], md["d"], color=PATINA, lw=1.8)
i = int(np.argmax(md["d"]))
b.plot(md["m"][i], md["d"][i], "o", color=PATINA, ms=6)
b.annotate(f"max $d$ = {md['d'][i]:.2f}", (md["m"][i], md["d"][i]),
           xytext=(6, 2), textcoords="offset points", color=INK, fontsize=8)
b.set_xlabel("edge count $m$  (density fixed by construction)")
b.set_ylabel(r"$|$Cohen $d|$, $\chi$ at matched density")
b.set_title("(b)  connectance-indexed: structure separates", loc="left")
save(fig, "fig1_confound")

# ---- fig 2: design guidance ------------------------------------------------
# (a) results/06_power_pegasus.txt: 499 perms, 50 realizations per cell
P = [24, 40, 80, 160]
power = [0.00, 0.08, 0.18, 0.66]
fpr = [0.08, 0.02, 0.02, 0.06]
# (b) results/16_multiyear_full.txt: 99 perms, 50 realizations per cell
T = [1, 3, 5, 10]
stacked = [0.26, 0.32, 0.44, 0.36]
averaged = [0.28, 0.44, 0.18, 0.10]

fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.5), sharey=True)
a.plot(P, power, "-o", color=COPPER, lw=1.8, ms=6)
a.plot(P, fpr, ":s", color=GRAY, lw=1.4, ms=5, markerfacecolor="white")
a.axhline(0.05, color=GRAY, lw=0.8, ls="--")
a.annotate("power (sharp threshold)", (P[3], power[3]), xytext=(-10, 0),
           textcoords="offset points", color=INK, fontsize=8, ha="right")
a.annotate("false-positive rate", (P[1], fpr[1]), xytext=(0, 7),
           textcoords="offset points", color=INK, fontsize=8)
a.annotate(r"$\alpha=.05$", (135, 0.05), xytext=(0, -11),
           textcoords="offset points", color=GRAY, fontsize=7.5)
a.set_xticks(P)
a.set_xlabel("plots $P$ (single survey)")
a.set_ylabel(r"rejection rate at $\alpha = .05$")
a.set_ylim(-0.03, 1.0)
a.set_title("(a)  power is bought with plots …", loc="left")

b.plot(T, stacked, "-o", color=COPPER, lw=1.8, ms=6)
b.plot(T, averaged, "--s", color=GRAY, lw=1.4, ms=5, markerfacecolor="white")
b.axhline(0.05, color=GRAY, lw=0.8, ls="--")
b.annotate("years stacked", (T[2], stacked[2]), xytext=(0, 7),
           textcoords="offset points", color=INK, fontsize=8)
b.annotate("years averaged", (T[2], averaged[2]), xytext=(0, -13),
           textcoords="offset points", color=INK, fontsize=8)
b.set_xticks(T)
b.set_xlabel("annual surveys $T$  (at $P=80$ plots)")
b.set_title("(b)  … only modestly with repeat surveys", loc="left")
save(fig, "fig2_design")

# ---- fig 3: the Everglades application -------------------------------------
import sys
sys.path.insert(0, os.path.join(ROOT, "experiments"))
glades = __import__("09_glades")          # recomputes nothing at import
C, zc = glades.curves(glades.env, "chi")  # 9 windows x 12 edge counts
Z = (C - C.mean(0)) / (C.std(0, ddof=1) + 1e-12)

fig, ax = plt.subplots(figsize=(6.6, 3.0))
for j in range(Z.shape[1]):
    ax.plot(zc, Z[:, j], color=GRAY, lw=0.7, alpha=0.35)
ax.plot(zc, Z.mean(1), color=PATINA, lw=2.2)
ax.annotate(r"mean standardized $\chi$", (zc[-3], Z.mean(1)[-3]),
            xytext=(0, 8), textcoords="offset points", color=INK, fontsize=8)

ax.axvspan(13.2, 19.9, color=PATINA, alpha=0.10, lw=0)
ax.axvspan(21.7, 34.1, color=PATINA, alpha=0.05, lw=0)
ax.axvline(15.1, color=PATINA, lw=1.0, ls=":")
ax.axvline(32.5, color=PATINA, lw=1.0, ls=":")
ax.axvline(10.0, color=GRAY, lw=1.0, ls="--")
ax.axvline(12.5, color=COPPER, lw=1.8)
ymax = ax.get_ylim()[1]
for x, txt, c in ((10.0, "10 criterion", GRAY),
                  (12.5, r"$\hat z$ = 12.5", COPPER),
                  (15.1, "TITAN decliners", PATINA),
                  (32.5, "TITAN increasers", PATINA)):
    ax.annotate(txt, (x, ymax), xytext=(3, -2), textcoords="offset points",
                color=c, fontsize=7.5, rotation=90, va="top")
ax.set_xscale("log")
ax.set_xticks([5, 10, 20, 50, 100])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_xlabel(r"window position: mean total phosphorus (μg/L, log scale)")
ax.set_ylabel(r"standardized $\chi(m)$ per window")
save(fig, "fig3_glades")

# ---- fig 4: the operating range (results/18_limits.txt) ------------------
fig, ax = plt.subplots(1, 4, figsize=(7.2, 2.1))
a, b, c, d = ax
a.plot([0.3, 0.5, 0.7, 0.9], [0.17, 0.23, 0.57, 0.80], "-o", color=COPPER, lw=1.8, ms=5)
a.set_xlabel("effect size $a$"); a.set_ylabel(r"power at $\alpha=.05$")
a.set_title("(a)  reorganization", loc="left", fontsize=8.5)
b.plot([12, 24, 60, 120], [0.13, 0.47, 0.67, 0.53], "-o", color=COPPER, lw=1.8, ms=5)
b.set_xlabel("taxa $S$"); b.set_title("(b)  richness", loc="left", fontsize=8.5)
lab = ["Gaussian", "80%", "40%", "15%", "pres/abs"]
vals = [0.60, 0.37, 0.13, 0.07, 0.13]
c.bar(range(5), vals, color=[PATINA] + [COPPER] * 3 + [GRAY], width=.68)
c.set_xticks(range(5)); c.set_xticklabels(lab, rotation=45, ha="right", fontsize=7)
c.set_title("(c)  data type", loc="left", fontsize=8.5)
d.bar([0, 1], [0.15, 0.07], color=[COPPER, PATINA], width=.6)
d.axhline(0.05, color=GRAY, lw=0.9, ls="--")
d.set_xticks([0, 1]); d.set_xticklabels(["plot-level", "within-region"], rotation=45, ha="right", fontsize=7)
d.set_ylabel("false-positive rate"); d.set_ylim(0, 0.2)
d.set_title("(d)  permutation", loc="left", fontsize=8.5)
for p_ in (a, b, c):
    p_.set_ylim(0, 0.9)
save(fig, "fig4_operating_range")
print("all figures written to paper/figs/")
