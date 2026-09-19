#!/usr/bin/env bash
# Public inputs for the ECS-on-omics proof of concept.  Nothing here is
# redistributed; this only downloads from the original sources.
set -e
cd "$(dirname "$0")"

# Arabidopsis 1001 Genomes expression, 727 accessions (Kawakatsu et al.)
[ -f ath_tx.tsv.gz ] || curl -sL -o ath_tx.tsv.gz \
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE80744&format=file&file=GSE80744%5Fath1001%5Ftx%5Fnorm%5F2016%2D04%2D21%2DUQ%5FgNorm%5FnormCounts%5Fk4%2Etsv%2Egz"

# Flowering time at 16C, from the Shiu Lab multi-omics repo
[ -f FT_16C.txt ] || gh api repos/ShiuLab/2024_Ath_GP/contents/Datasets/FT_16C.txt \
  --jq '.content' | base64 -d > FT_16C.txt

# Switchgrass drought and recovery (Meyer et al. 2014), 119 samples.
# The first input with a genuine ordered environmental axis.
if [ ! -d gse57887 ]; then
  curl -sL -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE57nnn/GSE57887/suppl/GSE57887_RAW.tar"
  mkdir -p gse57887 && tar -xf GSE57887_RAW.tar -C gse57887
fi
[ -f gse57887_samples.txt ] || curl -s \
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57887&targ=gsm&form=text&view=brief" \
  > gse57887_samples.txt

# Everglades macroinvertebrate benchmark (Baker & King 2010), shipped with the
# TITAN2 R package: 126 sites x 164 taxa, surface-water total phosphorus.
if [ ! -f glades_taxa.csv ]; then
  Rscript -e 'if (!requireNamespace("TITAN2", quietly=TRUE)) install.packages("TITAN2", repos="https://cloud.r-project.org");
    library(TITAN2); data(glades.taxa); data(glades.env);
    write.csv(glades.taxa, "glades_taxa.csv", row.names=FALSE);
    write.csv(glades.env,  "glades_env.csv",  row.names=FALSE)'
fi

ls -lh
