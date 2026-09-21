# Paper draft

Build:

    make          # main.pdf      11pt preprint
    make esa      # main_esa.pdf  double-spaced, 12pt, continuously line-numbered
    make tables   # regenerate tab_program.tex and tab_temporal.tex from ../results

`main.tex` selects its document class from `\ifdefined\ESAMODE`, so both
builds come from one source.

## What is generated, not written

`tab_program.tex` and `tab_temporal.tex` are produced by
`../experiments/32_tables.py` directly from the files in `../results`. No
application number in the manuscript is transcribed by hand. Edit the
generator, not the tables.

Figures come from `../experiments/13_paper_figures.py` into `figs/`.

## Structure

- §1 Introduction — the blind spot in taxon-aggregative threshold detection.
- §2 Methods — the connectance-indexed Euler characteristic summary; the
  pre-registered windowed estimator; the interaction test built on the
  intersection Euler characteristic profile; the two-stage procedure; the
  curvature attribution.
- §3 Results — simulation (estimator selection, power, operating range,
  closure, the instrument comparison, summary choice, complementarity with
  TITAN) then field data.
- §4 The registered cross-dataset ordering test.
- §5 The Butte Hill protocol, committed to in advance of the data.
- §6–8 Discussion, Conclusions, Open Research.

## Pre-registration

Every field dataset was registered before it was analyzed. The rules, the
dated amendments, the deviation records and the reproducibility correction
are in `../preregistration/`. Results that did not detect are reported
alongside those that did.
