"""Sensitivity of the pooled statistic to the edge-count grid density.

chi(m) is a cumulative sum, hence smooth in m; the pre-registered grid uses
12 evenly spaced counts. This check reruns the estimator with 6, 12, and 24
points -- on simulated THRESH/NULL gradients at P=80 (30 realizations) and
on the glades observed estimate -- to show the choice is immaterial.
"""
import importlib
import numpy as np

sim = importlib.import_module("06_gradient_estimator")

print("simulation, P=80, 30 realizations, split_mean")
print(f"{'grid':>5} {'THRESH stat':>12} {'NULL stat':>10} {'bias':>7} {'RMSE':>6}")
for G in (6, 12, 24):
    sim.M_GRID = np.unique(np.linspace(40, 900, G).astype(int))
    W, stride = sim.design(80)
    st, zh, nl = [], [], []
    for i in range(30):
        C, zc = sim.window_curves(*sim.community(80, "THRESH",
                                  np.random.default_rng(600 + i)), W, stride)
        z_, s_ = sim.est_split_mean(C, zc)
        st.append(s_); zh.append(z_)
        C, zc = sim.window_curves(*sim.community(80, "NULL",
                                  np.random.default_rng(600 + i)), W, stride)
        nl.append(sim.est_split_mean(C, zc)[1])
    zh = np.array(zh)
    print(f"{G:5d} {np.mean(st):12.2f} {np.mean(nl):10.2f} "
          f"{zh.mean() - sim.ZSTAR:+7.3f} {np.sqrt(((zh - sim.ZSTAR)**2).mean()):6.3f}")

gl = importlib.import_module("09_glades")
print("\nglades observed estimate")
print(f"{'grid':>5} {'z-hat':>7} {'stat':>6}")
for G in (6, 12, 24):
    S = gl.S
    gl.M_GRID = np.unique(np.linspace(int(np.ceil(2 * S / 3)),
                                      (S * (S - 1) // 2) // 2, G).astype(int))
    C, zc = gl.curves(gl.env, "chi")
    z_, s_ = gl.sim.est_split_mean(C, zc)
    print(f"{G:5d} {z_:7.1f} {s_:6.2f}")
