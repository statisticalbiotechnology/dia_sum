#!/bin/bash

echo "PyProhet export..."

pyprophet export --in=merged.osw --out MERGED.tsv --format legacy_merged --max_transition_pep 1.0 --ipf disable --max_global_peptide_qvalue 1.0 --max_rs_peakgroup_qvalue 1.0 --max_global_protein_qvalue 1.0

#pyprophet export --in=merged.scored.osw --out=legacy.tsv


