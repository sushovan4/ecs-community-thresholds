"""Generate the paper's application tables straight from the result files.

No number in Tables 5 and 6 is transcribed by hand: this script reads
results/28_mixup_*.txt (the interaction test, amendment 5), the matching
_reg.txt files (the registered 0.30-0.70 split grid), results/19_multi_*.txt
and results/22_temporal_*.txt (the curve statistic, first generation), and
writes paper/tab_program.tex and paper/tab_temporal.tex.

  python3 32_tables.py
"""
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")
P = os.path.join(ROOT, "paper")

# label, gradient, and which first-generation file holds the curve statistic
SPATIAL = [
    ("ELT", "Everglades LT insects", "days since dry", "17_everglades_lt_primary"),
    ("D1", "EPA rivers \\& streams", "conductivity", "19_multi_D1"),
    ("D2", "EPA streams, Appalachia", "conductivity", "19_multi_D2"),
    ("D3", "EMAP mine-drainage streams", "conductivity", "19_multi_D3"),
    ("D4", "EPA estuaries", "sediment copper", "19_multi_D4"),
    ("D5", "EPA lakes", "total phosphorus", "19_multi_D5"),
    ("D6", "Cedar Creek (N addition)", "N dose", "19_multi_D6"),
    ("D7", "Niwot alpine plants", "snow depth", "19_multi_D7"),
    ("D8", "CRMS Louisiana marsh", "station elevation", "19_multi_D8"),
]

TEMPORAL = [
    ("E2", "Everglades MWD SRS"), ("E3", "Everglades MWD TSL"),
    ("E4", "Everglades MWD WCA"), ("E6", "Everglades CERP SRS"),
    ("E8", "Everglades CERP WCA"), ("C1", "Cedar Creek field A"),
    ("C2", "Cedar Creek field B"), ("C3", "Cedar Creek field C"),
    ("C4", "Cedar Creek field D"), ("N1", "Niwot Saddle grid"),
    ("T2", "Everglades SRS, by occasion"),
    ("T3", "Everglades TSL, by occasion"),
    ("T4", "Everglades WCA, by occasion"),
]


def read(name):
    f = os.path.join(R, name if name.endswith(".txt") else name + ".txt")
    return open(f).read() if os.path.exists(f) else ""


def num(text, pattern, group=1):
    m = re.search(pattern, text)
    return float(m.group(group)) if m else None


FLOOR = 100          # units below which the interaction test over-rejects


def fmt(x, prec=3, bold_at=None, below_floor=False):
    """Bold marks a detection.  A p-value from a design below the floor of
    Section 3.7 is not evidence, so it is never bolded and is daggered."""
    if x is None:
        return "---"
    s = f"{x:.{prec}f}" if prec else f"{x:.4g}"
    if below_floor:
        return s + "$^\\dagger$"
    return f"\\textbf{{{s}}}" if bold_at is not None and x <= bold_at else s


def g(x):
    return "---" if x is None else f"{x:.4g}"


def spatial_rows():
    out = []
    for key, label, grad, curvefile in SPATIAL:
        mix, reg = read(f"28_mixup_{key}"), read(f"28_mixup_{key}_reg")
        cur = read(curvefile)
        units = num(mix, r"(\d+) units")
        spars = num(mix, r"sparsity ([0-9.]+)")
        p_mix = num(mix, r"p = ([0-9.]+)")
        p_rgn = num(mix, r"within-region \([^)]*\): p = ([0-9.]+)")
        p_grd = num(reg, r"p = ([0-9.]+)")
        z_cur = num(cur, r"z-hat = ([-0-9.e+]+)")
        p_cur = num(cur, r"z-hat = [-0-9.e+]+ ?d?\s+stat \S+\s+p = ([0-9.]+)")
        out.append((label, grad, units, spars, p_mix, p_rgn, p_grd, z_cur, p_cur))
    return out


def write_spatial():
    L = [r"\begin{table}[t]", r"\centering",
         r"\caption{The two-stage instrument on nine public gradients, all "
         r"pre-registered. $p$ is the interaction test of "
         r"Section~\ref{sec:mixup} over 199 unit permutations; "
         r"$p_{\text{reg}}$ permutes within regions, the test of record "
         r"where the gradient is spatially structured; $p_{\text{grid}}$ "
         r"repeats the test on the narrower split grid registered in "
         r"amendment~5. $\hat z$ and $p_\chi$ are the first-generation "
         r"curve statistic on the same matrix, reported unchanged. "
         r"``zeros'' is the fraction of zero entries, the pre-flight "
         r"diagnostic of Section~\ref{sec:limits}. A dagger marks a design "
         r"below the hundred-unit floor of Section~\ref{sec:mixloc}; its "
         r"$p$-values are reported but not counted as evidence.}",
         r"\label{tab:program}", r"\small",
         r"\begin{tabular}{llrrrrrrr}", r"\toprule",
         r"Dataset & Gradient & Units & zeros & $p$ & $p_{\text{reg}}$ "
         r"& $p_{\text{grid}}$ & $\hat z$ & $p_\chi$ \\", r"\midrule"]
    for (label, grad, units, spars, p_mix, p_rgn, p_grd, z_cur, p_cur) in spatial_rows():
        low = units is not None and units < FLOOR
        L.append(f"{label} & {grad} & "
                 f"{'---' if units is None else int(units)} & "
                 f"{'---' if spars is None else f'{spars:.2f}'} & "
                 f"{fmt(p_mix, 3, 0.05, low)} & {fmt(p_rgn, 3, 0.05, low)} & "
                 f"{fmt(p_grd, 3, 0.05, low)} & {g(z_cur)} & {fmt(p_cur, 3)} \\\\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(P, "tab_program.tex"), "w").write("\n".join(L) + "\n")
    return L


def write_temporal():
    L = [r"\begin{table}[t]", r"\centering",
         r"\caption{Time as the gradient: the interaction test on the "
         r"thirteen registered long-term series. Annual series take the "
         r"survey year as the unit; the three \texttt{T} series take the "
         r"survey occasion, five per year. $p_\chi$ is the "
         r"first-generation curve statistic on the same series, which "
         r"rejects only in T2 ($p_\chi = 0.045$). A dagger marks a series "
         r"below the hundred-unit design floor of "
         r"Section~\ref{sec:mixloc}, where the interaction test "
         r"over-rejects; those $p$-values are reported but are not "
         r"evidence, and are never bolded.}",
         r"\label{tab:temporal}", r"\small",
         r"\begin{tabular}{llrrrr}", r"\toprule",
         r"Series & Programme & Units & zeros & $p$ & $p_\chi$ \\",
         r"\midrule"]
    for key, label in TEMPORAL:
        mix, cur = read(f"28_mixup_{key}"), read(f"22_temporal_{key}")
        units = num(mix, r"(\d+) units")
        spars = num(mix, r"sparsity ([0-9.]+)")
        p_mix = num(mix, r"p = ([0-9.]+)")
        p_cur = num(cur, r"structural: .*p = ([0-9.]+)")
        low = units is not None and units < FLOOR
        L.append(f"{key} & {label} & "
                 f"{'---' if units is None else int(units)} & "
                 f"{'---' if spars is None else f'{spars:.2f}'} & "
                 f"{fmt(p_mix, 3, 0.05, low)} & {fmt(p_cur, 3)} \\\\")
    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    open(os.path.join(P, "tab_temporal.tex"), "w").write("\n".join(L) + "\n")
    return L


if __name__ == "__main__":
    for line in write_spatial() + [""] + write_temporal():
        print(line)
