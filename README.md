# PhyloChronoRank (PCR)

`PhyloChronoRank (PCR)` has two entry tabs.
It "amplifies" the signal already present in finished time trees by scoring them under a common set of diagnostics.

## [1_PCR_CUSTOM_DATING_PIPELINE_FROM_PHYLOGRAMS](1_PCR_CUSTOM_DATING_PIPELINE_FROM_PHYLOGRAMS/README.md)

Start from an unconstrained phylogram plus a calibration source. This tab runs `chronos`, `treePL`, and `RelTime` under one shared calibration-resolution step and writes a PCR-ready candidate set.

## [2_PCR_POSTFIT_METRICS](2_PCR_POSTFIT_METRICS/README.md)

Start from finished chronograms. This tab explains how PCR scores competing time trees with pulse, gap, rate, and optional uncertainty summaries, and it shows the bundled empirical examples.
