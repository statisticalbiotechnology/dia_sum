#!/bin/bash

echo "PyProhet merge..."

pyprophet merge --template=ecolihumanyeast_concat_mayu_IRR_cons_openswath_64w_fixed_curated_target_decoy.pqp --out=merged.osw *.osw

