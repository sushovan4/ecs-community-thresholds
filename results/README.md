# Result logs

Every number in the paper is read from a file here. Settings (permutations,
bootstrap resamples, realizations) are in each file's header.

| File | Paper location |
|---|---|
| `figures.json` | Fig. 1 (from experiments 01, 02) |
| `06_estimator_recovery_pegasus.txt` | Table 1 (100 realizations per design) |
| `06_power_pegasus.txt` | Sec. 3.4 power table, Fig. 2a |
| `06_bootstrap_pegasus.txt` | Sec. 3.4 sharpness calibration |
| `11_design_extras.txt` | Sec. 3.4 coverage; first (15-realization) repeat-survey runs |
| `16_multiyear_full.txt` | Sec. 3.4 repeat-survey power, Fig. 2b (50 realizations) |
| `08_summary_recovery_pegasus.txt`, `08_summary_rejection_pegasus.txt` | Table 2 |
| `07_titan.txt` | Sec. 3.6, Table 3 |
| `10_edivisive.txt`, `12_permanova.txt` | Sec. 3.5 baselines |
| `09_glades.txt`, `09_glades_titan.txt` | Sec. 4, Fig. 3 |
| `15_glades_anchored.txt` | Sec. 4, exploratory anchored test |
| `14_grid_sensitivity.txt` | Sec. 2.3 |
| `03_arabidopsis_199perm.txt`, `04_switchgrass_199perm.txt` | Sec. 4, omics stress tests |

Files without the `_pegasus` suffix in the 06/08 series are earlier,
lower-strength runs kept for provenance.

## Added 2026-09-19

| File | What it holds |
|---|---|
| `17_everglades_lt_primary.txt` | Everglades long-term panel: primary, anchored, within-region, robustness panel |
| `18_limits.txt` | the operating range: effect size, taxa, data type, threshold position, spatial confounding |
| `19_multi_D*.txt` | the frozen estimator on each public gradient (D1-D8) |
| `20_persistence_swap.txt` | persistence landscapes and total persistence against chi |
| `21_curvature_attribution.txt` | per-taxon curvature shares of the Everglades detection |
| `22_temporal_*.txt` | time as the gradient: ten annual series, three occasion-level series |
| `23_statistic_upgrade.txt` | pooling statistics compared (mean, Mahalanobis, CUSUM) |
| `24_window_scaling.txt` | window width against survey size |
| `25_titan_*.txt` | TITAN change points for the same datasets |
| `26_mixup_split.txt` | intersection-ECP split test against the registered statistic |
