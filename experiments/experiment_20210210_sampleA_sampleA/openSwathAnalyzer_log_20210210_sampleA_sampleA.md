(base) ptruong@planck:/hdd_14T/data/PXD002952/experiment_20210210_sampleA_sampleA$ ./experiment_20210210_sampleA_sampleA.sh 
Convert spectral library from .tsv to .TraML
Generate decoy
OpenSwathAnalyzer for .tsv output
Since neither rt_norm nor tr_irt is set, OpenSWATH will not use RT-transformation (rather a null transformation will be applied)
Progress of 'Load TraML file':
-- done [took 24.88 s (CPU), 24.88 s (Wall)] -- 
Loaded 13842 proteins, 91992 compounds with 551952 transitions.
Loading mzML file HYE110_TTOF6600_32var_lgillet_I160309_001_Pedro_Sample_A.mzML using readoptions cache
Progress of 'Loading metadata file HYE110_TTOF6600_32var_lgillet_I160309_001_Pedro_Sample_A.mzML':
Will analyze the metadata first to determine the number of SWATH windows and the window sizes.
Determined there to be 32 SWATH windows and in total 2197 MS1 spectra
Determined there to be 32 SWATH windows and in total 2197 MS1 spectra
-- done [took 01:14 m (CPU), 02:40 m (Wall)] -- 
Progress of 'Loading data file HYE110_TTOF6600_32var_lgillet_I160309_001_Pedro_Sample_A.mzML':
Read chromatogram while reading SWATH files, did not expect that!

  -- done [took 20:38 m (CPU), 17:43 m (Wall)] -- 
Will analyze 551952 transitions in total.

  Progress of 'Extracting and scoring transitions':
Use non-nested loop with 6 threads.
    3.03 %               Thread 0_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 0 out of 3)
Thread 3_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 0 out of 3)
Thread 1_0 will analyze 1828 compounds and 10968 transitions from SWATH 1 (batch 0 out of 1)
Thread 0_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 1 out of 3)
Thread 3_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 1 out of 3)
Thread 5_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 0 out of 3)
Thread 4_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 0 out of 3)
Thread 2_0 will analyze 2421 compounds and 14526 transitions from SWATH 2 (batch 0 out of 2)
Thread 1_0 will analyze 1828 compounds and 10968 transitions from SWATH 1 (batch 1 out of 1)
Thread 0_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 2 out of 3)
Thread 2_0 will analyze 2421 compounds and 14526 transitions from SWATH 2 (batch 1 out of 2)
  6.06 %               Thread 5_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 1 out of 3)
Thread 3_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 2 out of 3)
Thread 4_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 1 out of 3)
Thread 3_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 3 out of 3)
Thread 5_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 2 out of 3)
Thread 2_0 will analyze 2421 compounds and 14526 transitions from SWATH 2 (batch 2 out of 2)
  9.09 %               Thread 0_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 3 out of 3)
Thread 4_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 2 out of 3)
  12.12 %               Thread 5_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 3 out of 3)
Thread 4_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 3 out of 3)
  18.18 %               Thread 1_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 0 out of 3)
  21.21 %               Thread 1_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 1 out of 3)
Thread 1_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 2 out of 3)
Thread 1_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 3 out of 3)
  24.24 %               Thread 0_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 0 out of 4)
Thread 4_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 0 out of 3)
Thread 3_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 0 out of 3)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 0 out of 4)
Thread 0_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 1 out of 4)
Thread 3_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 1 out of 3)
Thread 5_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 0 out of 3)
Thread 4_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 1 out of 3)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 1 out of 4)
Thread 0_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 2 out of 4)
Thread 3_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 2 out of 3)
Thread 1_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 0 out of 3)
Thread 5_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 1 out of 3)
Thread 4_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 2 out of 3)
Thread 0_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 3 out of 4)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 2 out of 4)
Thread 3_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 3 out of 3)
Thread 1_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 1 out of 3)
Thread 5_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 2 out of 3)
Thread 4_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 3 out of 3)
              27.27 %               Thread 0_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 4 out of 4)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 3 out of 4)
Thread 1_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 2 out of 3)
  30.30 %               Thread 5_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 3 out of 3)
Thread 1_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 3 out of 3)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 4 out of 4)
      42.42 %               Thread 3_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 0 out of 3)
Thread 3_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 1 out of 3)
Thread 0_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 0 out of 3)
Thread 3_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 2 out of 3)
Thread 4_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 0 out of 3)
Thread 2_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 0 out of 3)
Thread 0_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 1 out of 3)
Thread 4_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 1 out of 3)
Thread 3_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 3 out of 3)
  45.45 %               Thread 2_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 1 out of 3)
Thread 4_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 2 out of 3)
Thread 0_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 2 out of 3)
Thread 2_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 2 out of 3)
Thread 4_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 3 out of 3)
Thread 0_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 3 out of 3)
                  51.52 %               Thread 2_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 3 out of 3)
                  54.55 %               Thread 5_0 will analyze 2765 compounds and 16590 transitions from SWATH 18 (batch 0 out of 2)
Thread 5_0 will analyze 2765 compounds and 16590 transitions from SWATH 18 (batch 1 out of 2)
Thread 5_0 will analyze 2765 compounds and 16590 transitions from SWATH 18 (batch 2 out of 2)
                  57.58 %               Thread 1_0 will analyze 2951 compounds and 17706 transitions from SWATH 19 (batch 0 out of 2)
Thread 1_0 will analyze 2951 compounds and 17706 transitions from SWATH 19 (batch 1 out of 2)
Thread 1_0 will analyze 2951 compounds and 17706 transitions from SWATH 19 (batch 2 out of 2)
                  60.61 %               Thread 3_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 0 out of 3)
Thread 3_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 1 out of 3)
Thread 3_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 2 out of 3)
Thread 4_0 will analyze 2888 compounds and 17328 transitions from SWATH 22 (batch 0 out of 2)
Thread 3_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 3 out of 3)
                63.64 %               Thread 0_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 0 out of 3)
Thread 2_0 will analyze 2756 compounds and 16536 transitions from SWATH 23 (batch 0 out of 2)
Thread 4_0 will analyze 2888 compounds and 17328 transitions from SWATH 22 (batch 1 out of 2)
Thread 0_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 1 out of 3)
Thread 2_0 will analyze 2756 compounds and 16536 transitions from SWATH 23 (batch 1 out of 2)
Thread 4_0 will analyze 2888 compounds and 17328 transitions from SWATH 22 (batch 2 out of 2)
Thread 2_0 will analyze 2756 compounds and 16536 transitions from SWATH 23 (batch 2 out of 2)
Thread 0_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 2 out of 3)
                                    69.70 %               Thread 0_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 3 out of 3)
                                    72.73 %               Thread 5_0 will analyze 2377 compounds and 14262 transitions from SWATH 24 (batch 0 out of 2)
Thread 5_0 will analyze 2377 compounds and 14262 transitions from SWATH 24 (batch 1 out of 2)
Thread 5_0 will analyze 2377 compounds and 14262 transitions from SWATH 24 (batch 2 out of 2)
                                    75.76 %               Thread 1_0 will analyze 2491 compounds and 14946 transitions from SWATH 25 (batch 0 out of 2)
Thread 1_0 will analyze 2491 compounds and 14946 transitions from SWATH 25 (batch 1 out of 2)
Thread 3_0 will analyze 2356 compounds and 14136 transitions from SWATH 26 (batch 0 out of 2)
Thread 1_0 will analyze 2491 compounds and 14946 transitions from SWATH 25 (batch 2 out of 2)
Thread 3_0 will analyze 2356 compounds and 14136 transitions from SWATH 26 (batch 1 out of 2)
                                      78.79 %               Thread 3_0 will analyze 2356 compounds and 14136 transitions from SWATH 26 (batch 2 out of 2)
                                      81.82 %               Thread 4_0 will analyze 2172 compounds and 13032 transitions from SWATH 27 (batch 0 out of 2)
Thread 2_0 will analyze 2029 compounds and 12174 transitions from SWATH 28 (batch 0 out of 2)
Thread 0_0 will analyze 1529 compounds and 9174 transitions from SWATH 29 (batch 0 out of 1)
Thread 4_0 will analyze 2172 compounds and 13032 transitions from SWATH 27 (batch 1 out of 2)
Thread 2_0 will analyze 2029 compounds and 12174 transitions from SWATH 28 (batch 1 out of 2)
Thread 4_0 will analyze 2172 compounds and 13032 transitions from SWATH 27 (batch 2 out of 2)
                                          84.85 %               Thread 0_0 will analyze 1529 compounds and 9174 transitions from SWATH 29 (batch 1 out of 1)
Thread 2_0 will analyze 2029 compounds and 12174 transitions from SWATH 28 (batch 2 out of 2)
                                          90.91 %               ./experiment_20210210_sampleA_sampleA.sh: line 13: 387400 Killed                  OpenSwathWorkflow -in HYE110_TTOF6600_32var_lgillet_I160309_001_Pedro_Sample_A.mzML -tr ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_var_curated_decoy.TraML -out_tsv osw_output.tsv -tempDirectory tmp -readOptions cacheWorkingInMemory -batchSize 1000 -Scoring:stop_report_after_feature -1 -min_upper_edge_dist 1 -extra_rt_extraction_window 100 -min_rsq 0.95 -min_coverage 0.6 -Scoring:Scores:use_dia_scores true -rt_extraction_window 600 -mz_extraction_window 30 -threads 6 -Scoring:TransitionGroupPicker:background_subtraction original
(base) ptruong@planck:/hdd_14T/data/PXD002952/experiment_20210210_sampleA_sampleA$ ls
ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_var_curated_decoy.TraML
ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_var_curated.TraML
ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_var_curated.tsv
experiment_20210210_sampleA_sampleA.sh
HYE110_TTOF6600_32var_lgillet_I160309_001_Pedro_Sample_A.mzML
HYE110_TTOF6600_32var_lgillet_I160309_003_Pedro_Sample_A.mzML
osw_output.tsv
tmp
(base) ptruong@planck:/hdd_14T/data/PXD002952/experiment_20210210_sampleA_sampleA$ vi experiment_20210210_sampleA_sampleA.sh 
(base) ptruong@planck:/hdd_14T/data/PXD002952/experiment_20210210_sampleA_sampleA$ ls
ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_var_curated_decoy.TraML
ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_var_curated.TraML
ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_var_curated.tsv
experiment_20210210_sampleA_sampleA.sh
HYE110_TTOF6600_32var_lgillet_I160309_001_Pedro_Sample_A.mzML
HYE110_TTOF6600_32var_lgillet_I160309_003_Pedro_Sample_A.mzML
osw_output.tsv
tmp
(base) ptruong@planck:/hdd_14T/data/PXD002952/experiment_20210210_sampleA_sampleA$ ./experiment_20210210_sampleA_sampleA.sh 
Convert spectral library from .tsv to .TraML
Generate decoy
OpenSwathAnalyzer for .tsv output
Since neither rt_norm nor tr_irt is set, OpenSWATH will not use RT-transformation (rather a null transformation will be applied)
Progress of 'Load TraML file':
-- done [took 25.30 s (CPU), 25.36 s (Wall)] -- 
Loaded 13842 proteins, 91992 compounds with 551952 transitions.
Loading mzML file HYE110_TTOF6600_32var_lgillet_I160309_001_Pedro_Sample_A.mzML using readoptions cache
Progress of 'Loading metadata file HYE110_TTOF6600_32var_lgillet_I160309_001_Pedro_Sample_A.mzML':
Will analyze the metadata first to determine the number of SWATH windows and the window sizes.
Determined there to be 32 SWATH windows and in total 2197 MS1 spectra
Determined there to be 32 SWATH windows and in total 2197 MS1 spectra
-- done [took 01:40 m (CPU), 02:45 m (Wall)] -- 
Progress of 'Loading data file HYE110_TTOF6600_32var_lgillet_I160309_001_Pedro_Sample_A.mzML':
Read chromatogram while reading SWATH files, did not expect that!

  -- done [took 19:27 m (CPU), 10:43 m (Wall)] -- 
Will analyze 551952 transitions in total.

  Progress of 'Extracting and scoring transitions':
Use non-nested loop with 6 threads.
    3.03 %               Thread 3_0 will analyze 2421 compounds and 14526 transitions from SWATH 2 (batch 0 out of 2)
Thread 2_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 0 out of 3)
Thread 0_0 will analyze 1828 compounds and 10968 transitions from SWATH 1 (batch 0 out of 1)
Thread 3_0 will analyze 2421 compounds and 14526 transitions from SWATH 2 (batch 1 out of 2)
Thread 1_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 0 out of 3)
Thread 5_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 0 out of 3)
Thread 4_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 0 out of 3)
Thread 2_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 1 out of 3)
Thread 0_0 will analyze 1828 compounds and 10968 transitions from SWATH 1 (batch 1 out of 1)
Thread 3_0 will analyze 2421 compounds and 14526 transitions from SWATH 2 (batch 2 out of 2)
6.06 %               Thread 1_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 1 out of 3)
Thread 2_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 2 out of 3)
Thread 4_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 1 out of 3)
  9.09 %               Thread 5_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 1 out of 3)
Thread 1_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 2 out of 3)
Thread 2_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 3 out of 3)
        12.12 %               Thread 4_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 2 out of 3)
Thread 5_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 2 out of 3)
Thread 1_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 3 out of 3)
                      15.15 %               Thread 4_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 3 out of 3)
Thread 5_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 3 out of 3)
                    21.21 %               Thread 3_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 0 out of 3)
Thread 3_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 1 out of 3)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 0 out of 4)
Thread 0_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 0 out of 3)
Thread 3_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 2 out of 3)
Thread 4_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 0 out of 3)
Thread 1_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 0 out of 4)
Thread 5_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 0 out of 3)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 1 out of 4)
Thread 3_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 3 out of 3)
Thread 4_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 1 out of 3)
  24.24 %               Thread 0_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 1 out of 3)
Thread 1_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 1 out of 4)
Thread 4_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 2 out of 3)
Thread 5_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 1 out of 3)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 2 out of 4)
Thread 0_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 2 out of 3)
Thread 1_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 2 out of 4)
Thread 4_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 3 out of 3)
Thread 3_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 0 out of 3)
Thread 5_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 2 out of 3)
27.27 %               Thread 1_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 3 out of 4)
Thread 0_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 3 out of 3)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 3 out of 4)
Thread 3_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 1 out of 3)
Thread 1_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 4 out of 4)
  30.30 %               Thread 5_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 3 out of 3)
Thread 2_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 4 out of 4)
Thread 3_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 2 out of 3)
          39.39 %               Thread 3_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 3 out of 3)
          42.42 %               Thread 1_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 0 out of 3)
Thread 5_0 will analyze 2765 compounds and 16590 transitions from SWATH 18 (batch 0 out of 2)
Thread 4_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 0 out of 3)
Thread 1_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 1 out of 3)
Thread 0_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 0 out of 3)
Thread 2_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 0 out of 3)
Thread 4_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 1 out of 3)
Thread 5_0 will analyze 2765 compounds and 16590 transitions from SWATH 18 (batch 1 out of 2)
Thread 3_0 will analyze 2951 compounds and 17706 transitions from SWATH 19 (batch 0 out of 2)
Thread 0_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 1 out of 3)
Thread 1_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 2 out of 3)
Thread 5_0 will analyze 2765 compounds and 16590 transitions from SWATH 18 (batch 2 out of 2)
Thread 4_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 2 out of 3)
Thread 2_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 1 out of 3)
    45.45 %               Thread 3_0 will analyze 2951 compounds and 17706 transitions from SWATH 19 (batch 1 out of 2)
Thread 0_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 2 out of 3)
Thread 4_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 3 out of 3)
      48.48 %               Thread 1_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 3 out of 3)
Thread 2_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 2 out of 3)
Thread 3_0 will analyze 2951 compounds and 17706 transitions from SWATH 19 (batch 2 out of 2)
Thread 0_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 3 out of 3)
Thread 2_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 3 out of 3)
                        60.61 %               Thread 5_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 0 out of 3)
Thread 1_0 will analyze 2888 compounds and 17328 transitions from SWATH 22 (batch 0 out of 2)
Thread 4_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 0 out of 3)
Thread 5_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 1 out of 3)
Thread 3_0 will analyze 2756 compounds and 16536 transitions from SWATH 23 (batch 0 out of 2)
Thread 1_0 will analyze 2888 compounds and 17328 transitions from SWATH 22 (batch 1 out of 2)
Thread 2_0 will analyze 2491 compounds and 14946 transitions from SWATH 25 (batch 0 out of 2)
Thread 4_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 1 out of 3)
Thread 5_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 2 out of 3)
Thread 0_0 will analyze 2377 compounds and 14262 transitions from SWATH 24 (batch 0 out of 2)
Thread 3_0 will analyze 2756 compounds and 16536 transitions from SWATH 23 (batch 1 out of 2)
Thread 1_0 will analyze 2888 compounds and 17328 transitions from SWATH 22 (batch 2 out of 2)
Thread 4_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 2 out of 3)
                        63.64 %               Thread 3_0 will analyze 2756 compounds and 16536 transitions from SWATH 23 (batch 2 out of 2)
Thread 2_0 will analyze 2491 compounds and 14946 transitions from SWATH 25 (batch 1 out of 2)
Thread 5_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 3 out of 3)
          66.67 %               Thread 0_0 will analyze 2377 compounds and 14262 transitions from SWATH 24 (batch 1 out of 2)
              69.70 %               Thread 4_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 3 out of 3)
Thread 2_0 will analyze 2491 compounds and 14946 transitions from SWATH 25 (batch 2 out of 2)
              75.76 %               Thread 0_0 will analyze 2377 compounds and 14262 transitions from SWATH 24 (batch 2 out of 2)
              78.79 %               Thread 1_0 will analyze 2356 compounds and 14136 transitions from SWATH 26 (batch 0 out of 2)
Thread 5_0 will analyze 2172 compounds and 13032 transitions from SWATH 27 (batch 0 out of 2)
Thread 4_0 will analyze 1529 compounds and 9174 transitions from SWATH 29 (batch 0 out of 1)
Thread 1_0 will analyze 2356 compounds and 14136 transitions from SWATH 26 (batch 1 out of 2)
Thread 3_0 will analyze 2029 compounds and 12174 transitions from SWATH 28 (batch 0 out of 2)
Thread 2_0 will analyze 1577 compounds and 9462 transitions from SWATH 30 (batch 0 out of 1)
Thread 0_0 will analyze 1201 compounds and 7206 transitions from SWATH 31 (batch 0 out of 1)
Thread 1_0 will analyze 2356 compounds and 14136 transitions from SWATH 26 (batch 2 out of 2)
Thread 4_0 will analyze 1529 compounds and 9174 transitions from SWATH 29 (batch 1 out of 1)
Thread 5_0 will analyze 2172 compounds and 13032 transitions from SWATH 27 (batch 1 out of 2)
Thread 3_0 will analyze 2029 compounds and 12174 transitions from SWATH 28 (batch 1 out of 2)
      84.85 %               Thread 2_0 will analyze 1577 compounds and 9462 transitions from SWATH 30 (batch 1 out of 1)
Thread 3_0 will analyze 2029 compounds and 12174 transitions from SWATH 28 (batch 2 out of 2)
Thread 5_0 will analyze 2172 compounds and 13032 transitions from SWATH 27 (batch 2 out of 2)
  90.91 %               Thread 0_0 will analyze 1201 compounds and 7206 transitions from SWATH 31 (batch 1 out of 1)
96.97 %               Thread 1_0 will analyze 546 compounds and 3276 transitions from SWATH 32 (batch 0 out of 1)
Thread 1_0 will analyze 546 compounds and 3276 transitions from SWATH 32 (batch 1 out of 1)
-- done [took 08:54 m (CPU), 06:13 m (Wall)] -- 
OpenSwathWorkflow took 20:08 m (wall), 30:29 m (CPU), 02:09 m (system), 28:19 m (user); Peak Memory Usage: 22366 MB.
Since neither rt_norm nor tr_irt is set, OpenSWATH will not use RT-transformation (rather a null transformation will be applied)
Progress of 'Load TraML file':
-- done [took 25.55 s (CPU), 26.13 s (Wall)] -- 
Loaded 13842 proteins, 91992 compounds with 551952 transitions.
Loading mzML file HYE110_TTOF6600_32var_lgillet_I160309_003_Pedro_Sample_A.mzML using readoptions cache
Progress of 'Loading metadata file HYE110_TTOF6600_32var_lgillet_I160309_003_Pedro_Sample_A.mzML':
Will analyze the metadata first to determine the number of SWATH windows and the window sizes.
Determined there to be 32 SWATH windows and in total 2197 MS1 spectra
Determined there to be 32 SWATH windows and in total 2197 MS1 spectra
-- done [took 01:42 m (CPU), 02:52 m (Wall)] -- 
Progress of 'Loading data file HYE110_TTOF6600_32var_lgillet_I160309_003_Pedro_Sample_A.mzML':
Read chromatogram while reading SWATH files, did not expect that!

  -- done [took 20:49 m (CPU), 18:08 m (Wall)] -- 
Will analyze 551952 transitions in total.

  Progress of 'Extracting and scoring transitions':
Use non-nested loop with 6 threads.
    3.03 %               Thread 2_0 will analyze 1828 compounds and 10968 transitions from SWATH 1 (batch 0 out of 1)
Thread 2_0 will analyze 1828 compounds and 10968 transitions from SWATH 1 (batch 1 out of 1)
Thread 5_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 0 out of 3)
    6.06 %               Thread 3_0 will analyze 2421 compounds and 14526 transitions from SWATH 2 (batch 0 out of 2)
Thread 5_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 1 out of 3)
Thread 1_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 0 out of 3)
Thread 4_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 0 out of 3)
Thread 3_0 will analyze 2421 compounds and 14526 transitions from SWATH 2 (batch 1 out of 2)
Thread 5_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 2 out of 3)
Thread 0_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 0 out of 3)
Thread 4_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 1 out of 3)
Thread 1_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 1 out of 3)
Thread 3_0 will analyze 2421 compounds and 14526 transitions from SWATH 2 (batch 2 out of 2)
                    9.09 %               Thread 4_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 2 out of 3)
Thread 1_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 2 out of 3)
Thread 5_0 will analyze 3261 compounds and 19566 transitions from SWATH 3 (batch 3 out of 3)
Thread 2_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 0 out of 3)
                                            12.12 %               Thread 0_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 1 out of 3)
Thread 4_0 will analyze 3192 compounds and 19152 transitions from SWATH 4 (batch 3 out of 3)
Thread 1_0 will analyze 3945 compounds and 23670 transitions from SWATH 5 (batch 3 out of 3)
                                                                                          15.15 %               Thread 2_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 1 out of 3)
Thread 0_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 2 out of 3)
                                                                                                18.18 %               Thread 2_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 2 out of 3)
Thread 0_0 will analyze 3443 compounds and 20658 transitions from SWATH 6 (batch 3 out of 3)
Thread 2_0 will analyze 3476 compounds and 20856 transitions from SWATH 7 (batch 3 out of 3)
                                                                                                                                                                                                        24.24 %               Thread 3_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 0 out of 3)
Thread 5_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 0 out of 4)
Thread 3_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 1 out of 3)
Thread 2_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 0 out of 3)
Thread 0_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 0 out of 3)
Thread 1_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 0 out of 3)
Thread 5_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 1 out of 4)
Thread 4_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 0 out of 4)
Thread 2_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 1 out of 3)
Thread 3_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 2 out of 3)
Thread 5_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 2 out of 4)
Thread 0_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 1 out of 3)
Thread 2_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 2 out of 3)
Thread 1_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 1 out of 3)
Thread 4_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 1 out of 4)
Thread 5_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 3 out of 4)
Thread 2_0 will analyze 3615 compounds and 21690 transitions from SWATH 13 (batch 3 out of 3)
Thread 3_0 will analyze 3872 compounds and 23232 transitions from SWATH 8 (batch 3 out of 3)
Thread 5_0 will analyze 4088 compounds and 24528 transitions from SWATH 9 (batch 4 out of 4)
    27.27 %               Thread 0_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 2 out of 3)
Thread 1_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 2 out of 3)
Thread 4_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 2 out of 4)
    33.33 %               Thread 1_0 will analyze 3759 compounds and 22554 transitions from SWATH 11 (batch 3 out of 3)
Thread 0_0 will analyze 3636 compounds and 21816 transitions from SWATH 12 (batch 3 out of 3)
Thread 4_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 3 out of 4)
                              39.39 %               Thread 4_0 will analyze 4131 compounds and 24786 transitions from SWATH 10 (batch 4 out of 4)
                              42.42 %               Thread 5_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 0 out of 3)
Thread 5_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 1 out of 3)
Thread 5_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 2 out of 3)
Thread 2_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 0 out of 3)
Thread 1_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 0 out of 3)
Thread 5_0 will analyze 3104 compounds and 18624 transitions from SWATH 15 (batch 3 out of 3)
                              45.45 %               Thread 2_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 1 out of 3)
Thread 1_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 1 out of 3)
Thread 3_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 0 out of 3)
Thread 2_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 2 out of 3)
Thread 1_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 2 out of 3)
Thread 2_0 will analyze 3239 compounds and 19434 transitions from SWATH 14 (batch 3 out of 3)
Thread 3_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 1 out of 3)
                                                  48.48 %               Thread 1_0 will analyze 3319 compounds and 19914 transitions from SWATH 17 (batch 3 out of 3)
                                                          51.52 %               Thread 3_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 2 out of 3)
Thread 0_0 will analyze 2765 compounds and 16590 transitions from SWATH 18 (batch 0 out of 2)
Thread 3_0 will analyze 3623 compounds and 21738 transitions from SWATH 16 (batch 3 out of 3)
Thread 0_0 will analyze 2765 compounds and 16590 transitions from SWATH 18 (batch 1 out of 2)
                                                    54.55 %               Thread 4_0 will analyze 2951 compounds and 17706 transitions from SWATH 19 (batch 0 out of 2)
Thread 0_0 will analyze 2765 compounds and 16590 transitions from SWATH 18 (batch 2 out of 2)
Thread 4_0 will analyze 2951 compounds and 17706 transitions from SWATH 19 (batch 1 out of 2)
                                                        57.58 %               Thread 4_0 will analyze 2951 compounds and 17706 transitions from SWATH 19 (batch 2 out of 2)
                                                        60.61 %               Thread 1_0 will analyze 2888 compounds and 17328 transitions from SWATH 22 (batch 0 out of 2)
Thread 5_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 0 out of 3)
Thread 2_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 0 out of 3)
Thread 3_0 will analyze 2756 compounds and 16536 transitions from SWATH 23 (batch 0 out of 2)
Thread 5_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 1 out of 3)
Thread 1_0 will analyze 2888 compounds and 17328 transitions from SWATH 22 (batch 1 out of 2)
Thread 0_0 will analyze 2377 compounds and 14262 transitions from SWATH 24 (batch 0 out of 2)
Thread 2_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 1 out of 3)
Thread 4_0 will analyze 2491 compounds and 14946 transitions from SWATH 25 (batch 0 out of 2)
Thread 1_0 will analyze 2888 compounds and 17328 transitions from SWATH 22 (batch 2 out of 2)
Thread 2_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 2 out of 3)
Thread 3_0 will analyze 2756 compounds and 16536 transitions from SWATH 23 (batch 1 out of 2)
Thread 5_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 2 out of 3)
                            63.64 %               Thread 2_0 will analyze 3017 compounds and 18102 transitions from SWATH 21 (batch 3 out of 3)
Thread 0_0 will analyze 2377 compounds and 14262 transitions from SWATH 24 (batch 1 out of 2)
Thread 4_0 will analyze 2491 compounds and 14946 transitions from SWATH 25 (batch 1 out of 2)
Thread 3_0 will analyze 2756 compounds and 16536 transitions from SWATH 23 (batch 2 out of 2)
Thread 5_0 will analyze 3081 compounds and 18486 transitions from SWATH 20 (batch 3 out of 3)
                                    69.70 %               Thread 4_0 will analyze 2491 compounds and 14946 transitions from SWATH 25 (batch 2 out of 2)
                                            72.73 %               Thread 0_0 will analyze 2377 compounds and 14262 transitions from SWATH 24 (batch 2 out of 2)
                                        78.79 %               Thread 3_0 will analyze 1529 compounds and 9174 transitions from SWATH 29 (batch 0 out of 1)
Thread 1_0 will analyze 2356 compounds and 14136 transitions from SWATH 26 (batch 0 out of 2)
Thread 2_0 will analyze 2172 compounds and 13032 transitions from SWATH 27 (batch 0 out of 2)
Thread 5_0 will analyze 2029 compounds and 12174 transitions from SWATH 28 (batch 0 out of 2)
Thread 3_0 will analyze 1529 compounds and 9174 transitions from SWATH 29 (batch 1 out of 1)
Thread 1_0 will analyze 2356 compounds and 14136 transitions from SWATH 26 (batch 1 out of 2)
                                                  81.82 %               Thread 2_0 will analyze 2172 compounds and 13032 transitions from SWATH 27 (batch 1 out of 2)
Thread 5_0 will analyze 2029 compounds and 12174 transitions from SWATH 28 (batch 1 out of 2)
Thread 1_0 will analyze 2356 compounds and 14136 transitions from SWATH 26 (batch 2 out of 2)
                                                  84.85 %               Thread 2_0 will analyze 2172 compounds and 13032 transitions from SWATH 27 (batch 2 out of 2)
                                                    87.88 %               Thread 5_0 will analyze 2029 compounds and 12174 transitions from SWATH 28 (batch 2 out of 2)
Thread 4_0 will analyze 1577 compounds and 9462 transitions from SWATH 30 (batch 0 out of 1)
Thread 0_0 will analyze 1201 compounds and 7206 transitions from SWATH 31 (batch 0 out of 1)
Thread 4_0 will analyze 1577 compounds and 9462 transitions from SWATH 30 (batch 1 out of 1)
Thread 0_0 will analyze 1201 compounds and 7206 transitions from SWATH 31 (batch 1 out of 1)
                                                  96.97 %               Thread 3_0 will analyze 546 compounds and 3276 transitions from SWATH 32 (batch 0 out of 1)
Thread 3_0 will analyze 546 compounds and 3276 transitions from SWATH 32 (batch 1 out of 1)
                                                  100.00 %               
                                                -- done [took 08:29 m (CPU), 06:22 m (Wall)] -- 
OpenSwathWorkflow took 27:49 m (wall), 31:27 m (CPU), 02:24 m (system), 29:02 m (user); Peak Memory Usage: 22251 MB.
