# PhyloChronoRank (PCR)

`PhyloChronoRank (PCR)` includes tools to "amplify" the signal in phylograms into alternative chronograms, then score those chronograms under a common set of diagnostics. It has two entry tabs:

## [1_PCR_CUSTOM_DATING_PIPELINE_FROM_PHYLOGRAMS](1_PCR_CUSTOM_DATING_PIPELINE_FROM_PHYLOGRAMS/README.md)

Start from an unconstrained phylogram plus a calibration source. This tab runs `chronos`, `treePL`, and `RelTime` under one shared calibration-resolution step to produce a set of candidate chronograms for comparison in tab 2.

## [2_PCR_POSTFIT_METRICS](2_PCR_POSTFIT_METRICS/README.md)

Start from finished chronograms. This tab asks which candidate time tree is the most biologically defensible: which one best preserves branching structure, stays closest to the calibration evidence, and avoids implausible rate behavior. It also shows the bundled empirical examples.
