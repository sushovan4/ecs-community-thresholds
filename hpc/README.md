# HPC runs (GWU Pegasus)

`pegasus_full.sbatch` reproduces the paper's simulation results at publication
strength (NPERM=499, 50-100 realisations per cell) and adds the P=160 design
the local runs skip — the interesting open question there is whether the
bootstrap CI-width sharpness diagnostic, which failed at P=80
(results/06_estimator_slow.txt), comes alive with narrower windows.

Local reduced-strength results are committed under `results/`; a Pegasus rerun
should replace them via the same redirection into `results/`, noting the
settings in the log header. Seeds are fixed in the scripts, so same-setting
reruns are exact reproductions.
