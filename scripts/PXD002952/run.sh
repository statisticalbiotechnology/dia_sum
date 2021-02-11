#!/bin/bash

for filename in *.mzML
do
  OpenSwathWorkflow -in $filename -tr ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_fixed_curated_decoy.TraML -out_tsv osw_output_$filename.tsv -readOptions cacheWorkingInMemory -batchSize 1000 -min_upper_edge_dist 1 -extra_rt_extraction_window 100 -min_rsq 0.95 -min_coverage 0.6 -rt_extraction_window 600 -mz_extraction_window 30 -threads 6 -ppm

done;



#OpenSwathWorkflow -in ../test_HYE124_TTOF6600_32fix_lgillet_I150211_004_test_sample.mzML -tr ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_fixed_curated_decoy.TraML -out_tsv osw_output.tsv -readOptions cacheWorkingInMemory -batchSize 1000 -min_upper_edge_dist 1 -extra_rt_extraction_window 100 -min_rsq 0.95 -min_coverage 0.6 -rt_extraction_window 600 -mz_extraction_window 30 -threads 6 -ppm




