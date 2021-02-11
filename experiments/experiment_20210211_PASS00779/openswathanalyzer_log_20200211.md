
(base) ptruong@planck:/hdd_14T/data/PASS00779/ftp.peptideatlas.org/experiment_20200211$ ./openswathanalyzer_run.sh 
.mzMl OpenSwathAbakyzer
Validate provided Swath windows file:
Read Swath window header: 'start	end'
Read Swath window file with 32 SWATH windows.
Progress of 'Load TSV file':
-- done [took 5.17 s (CPU), 5.17 s (Wall)] -- 
Loaded 8029 proteins, 95466 compounds with 477330 transitions.
Loading mzML file olgas_K121026_001_SW_Wayne_R1_d00-Mtb_Wayne_01.mzML using readoptions cache
Progress of 'Loading metadata file olgas_K121026_001_SW_Wayne_R1_d00-Mtb_Wayne_01.mzML':
Will analyze the metadata first to determine the number of SWATH windows and the window sizes.
Determined there to be 32 SWATH windows and in total 2277 MS1 spectra
Determined there to be 32 SWATH windows and in total 2277 MS1 spectra
-- done [took 14.09 s (CPU), 14.11 s (Wall)] -- 
Progress of 'Loading data file olgas_K121026_001_SW_Wayne_R1_d00-Mtb_Wayne_01.mzML':
Read chromatogram while reading SWATH files, did not expect that!

  -- done [took 58.94 s (CPU), 59.20 s (Wall)] -- 
Read Swath window header: 'start	end'
Read Swath window file with 32 SWATH windows.
Re-annotate from file: SWATH 400 / 425 (raw data) is annotated via swath_windows_file with 400 / 424.5
Re-annotate from file: SWATH 424 / 450 (raw data) is annotated via swath_windows_file with 424.5 / 449.5
Re-annotate from file: SWATH 449 / 475 (raw data) is annotated via swath_windows_file with 449.5 / 474.5
Re-annotate from file: SWATH 474 / 500 (raw data) is annotated via swath_windows_file with 474.5 / 499.5
Re-annotate from file: SWATH 499 / 525 (raw data) is annotated via swath_windows_file with 499.5 / 524.5
Re-annotate from file: SWATH 524 / 550 (raw data) is annotated via swath_windows_file with 524.5 / 549.5
Re-annotate from file: SWATH 549 / 575 (raw data) is annotated via swath_windows_file with 549.5 / 574.5
Re-annotate from file: SWATH 574 / 600 (raw data) is annotated via swath_windows_file with 574.5 / 599.5
Re-annotate from file: SWATH 599 / 625 (raw data) is annotated via swath_windows_file with 599.5 / 624.5
Re-annotate from file: SWATH 624 / 650 (raw data) is annotated via swath_windows_file with 624.5 / 649.5
Re-annotate from file: SWATH 649 / 675 (raw data) is annotated via swath_windows_file with 649.5 / 674.5
Re-annotate from file: SWATH 674 / 700 (raw data) is annotated via swath_windows_file with 674.5 / 699.5
Re-annotate from file: SWATH 699 / 725 (raw data) is annotated via swath_windows_file with 699.5 / 724.5
Re-annotate from file: SWATH 724 / 750 (raw data) is annotated via swath_windows_file with 724.5 / 749.5
Re-annotate from file: SWATH 749 / 775 (raw data) is annotated via swath_windows_file with 749.5 / 774.5
Re-annotate from file: SWATH 774 / 800 (raw data) is annotated via swath_windows_file with 774.5 / 799.5
Re-annotate from file: SWATH 799 / 825 (raw data) is annotated via swath_windows_file with 799.5 / 824.5
Re-annotate from file: SWATH 824 / 850 (raw data) is annotated via swath_windows_file with 824.5 / 849.5
Re-annotate from file: SWATH 849 / 875 (raw data) is annotated via swath_windows_file with 849.5 / 874.5
Re-annotate from file: SWATH 874 / 900 (raw data) is annotated via swath_windows_file with 874.5 / 899.5
Re-annotate from file: SWATH 899 / 925 (raw data) is annotated via swath_windows_file with 899.5 / 924.5
Re-annotate from file: SWATH 924 / 950 (raw data) is annotated via swath_windows_file with 924.5 / 949.5
Re-annotate from file: SWATH 949 / 975 (raw data) is annotated via swath_windows_file with 949.5 / 974.5
Re-annotate from file: SWATH 974 / 1000 (raw data) is annotated via swath_windows_file with 974.5 / 999.5
Re-annotate from file: SWATH 999 / 1025 (raw data) is annotated via swath_windows_file with 999.5 / 1024.5
Re-annotate from file: SWATH 1024 / 1050 (raw data) is annotated via swath_windows_file with 1024.5 / 1049.5
Re-annotate from file: SWATH 1049 / 1075 (raw data) is annotated via swath_windows_file with 1049.5 / 1074.5
Re-annotate from file: SWATH 1074 / 1100 (raw data) is annotated via swath_windows_file with 1074.5 / 1099.5
Re-annotate from file: SWATH 1099 / 1125 (raw data) is annotated via swath_windows_file with 1099.5 / 1124.5
Re-annotate from file: SWATH 1124 / 1150 (raw data) is annotated via swath_windows_file with 1124.5 / 1149.5
Re-annotate from file: SWATH 1149 / 1175 (raw data) is annotated via swath_windows_file with 1149.5 / 1174.5
Re-annotate from file: SWATH 1174 / 1200 (raw data) is annotated via swath_windows_file with 1174.5 / 1200
Will load iRT transitions and try to find iRT peptides

  Progress of 'Load TraML file':

  -- done [took 0.04 s (CPU), 0.05 s (Wall)] -- 

  Progress of 'Extract iRT chromatograms':

  -- done [took 0.37 s (CPU), 0.36 s (Wall)] -- 

  Progress of 'Retention time normalization':
Will analyse 20 peptides with a total of 100 transitions 
WARNING in SignalToNoiseEstimatorMedian: 1.05402% of all Signal-to-Noise estimates are too high, because the median was found in the rightmost histogram-bin. You should consider increasing 'max_intensity' (and maybe 'bin_count' with it, to keep bin width reasonable)
rsq: 0.427824 points: 20
rsq: 0.923524 points: 19
rsq: 0.999565 points: 18

  -- done [took 0.65 s (CPU), 0.66 s (Wall)] -- 
Will analyze 477330 transitions in total.

  Progress of 'Extracting and scoring transitions':
Use non-nested loop with 1 threads.
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 0 out of 3)
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 1 out of 3)
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 2 out of 3)
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 3 out of 3)
    6.06 %               Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 0 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 1 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 2 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 3 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 4 out of 4)
    9.09 %               Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 0 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 1 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 2 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 3 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 4 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 5 out of 5)
    12.12 %               Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 0 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 1 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 2 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 3 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 4 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 5 out of 5)
    15.15 %               Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 0 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 1 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 2 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 3 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 4 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 5 out of 5)
    18.18 %               Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 0 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 1 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 2 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 3 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 4 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 5 out of 5)
    21.21 %               Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 0 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 1 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 2 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 3 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 4 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 5 out of 5)
    24.24 %               Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 0 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 1 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 2 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 3 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 4 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 5 out of 5)
    27.27 %               Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 0 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 1 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 2 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 3 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 4 out of 4)
    30.30 %               Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 0 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 1 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 2 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 3 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 4 out of 4)
    33.33 %               Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 0 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 1 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 2 out of 4)

Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 3 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 4 out of 4)
    36.36 %               Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 0 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 1 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 2 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 3 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 4 out of 4)
    39.39 %               Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 0 out of 3)
Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 1 out of 3)
Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 2 out of 3)
Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 3 out of 3)
    42.42 %               Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 0 out of 3)
Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 1 out of 3)
Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 2 out of 3)
Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 3 out of 3)
    45.45 %               Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 0 out of 3)
Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 1 out of 3)
Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 2 out of 3)
Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 3 out of 3)
    48.48 %               Thread 0_0 will analyze 2972 compounds and 14860 transitions from SWATH 16 (batch 0 out of 2)
Thread 0_0 will analyze 2972 compounds and 14860 transitions from SWATH 16 (batch 1 out of 2)
Thread 0_0 will analyze 2972 compounds and 14860 transitions from SWATH 16 (batch 2 out of 2)
    51.52 %               Thread 0_0 will analyze 2689 compounds and 13445 transitions from SWATH 17 (batch 0 out of 2)
Thread 0_0 will analyze 2689 compounds and 13445 transitions from SWATH 17 (batch 1 out of 2)
Thread 0_0 will analyze 2689 compounds and 13445 transitions from SWATH 17 (batch 2 out of 2)
    54.55 %               Thread 0_0 will analyze 2391 compounds and 11955 transitions from SWATH 18 (batch 0 out of 2)
Thread 0_0 will analyze 2391 compounds and 11955 transitions from SWATH 18 (batch 1 out of 2)
Thread 0_0 will analyze 2391 compounds and 11955 transitions from SWATH 18 (batch 2 out of 2)
    57.58 %               Thread 0_0 will analyze 2171 compounds and 10855 transitions from SWATH 19 (batch 0 out of 2)
Thread 0_0 will analyze 2171 compounds and 10855 transitions from SWATH 19 (batch 1 out of 2)
Thread 0_0 will analyze 2171 compounds and 10855 transitions from SWATH 19 (batch 2 out of 2)
    60.61 %               Thread 0_0 will analyze 2019 compounds and 10095 transitions from SWATH 20 (batch 0 out of 2)
Thread 0_0 will analyze 2019 compounds and 10095 transitions from SWATH 20 (batch 1 out of 2)
Thread 0_0 will analyze 2019 compounds and 10095 transitions from SWATH 20 (batch 2 out of 2)
    63.64 %               Thread 0_0 will analyze 1829 compounds and 9145 transitions from SWATH 21 (batch 0 out of 1)
Thread 0_0 will analyze 1829 compounds and 9145 transitions from SWATH 21 (batch 1 out of 1)
    66.67 %               Thread 0_0 will analyze 1628 compounds and 8140 transitions from SWATH 22 (batch 0 out of 1)
Thread 0_0 will analyze 1628 compounds and 8140 transitions from SWATH 22 (batch 1 out of 1)
    69.70 %               Thread 0_0 will analyze 1526 compounds and 7630 transitions from SWATH 23 (batch 0 out of 1)
Thread 0_0 will analyze 1526 compounds and 7630 transitions from SWATH 23 (batch 1 out of 1)
    72.73 %               Thread 0_0 will analyze 1371 compounds and 6855 transitions from SWATH 24 (batch 0 out of 1)
Thread 0_0 will analyze 1371 compounds and 6855 transitions from SWATH 24 (batch 1 out of 1)
    75.76 %               Thread 0_0 will analyze 1185 compounds and 5925 transitions from SWATH 25 (batch 0 out of 1)
Thread 0_0 will analyze 1185 compounds and 5925 transitions from SWATH 25 (batch 1 out of 1)
    78.79 %               Thread 0_0 will analyze 1106 compounds and 5530 transitions from SWATH 26 (batch 0 out of 1)
Thread 0_0 will analyze 1106 compounds and 5530 transitions from SWATH 26 (batch 1 out of 1)
    81.82 %               Thread 0_0 will analyze 889 compounds and 4445 transitions from SWATH 27 (batch 0 out of 1)
Thread 0_0 will analyze 889 compounds and 4445 transitions from SWATH 27 (batch 1 out of 1)
    84.85 %               Thread 0_0 will analyze 766 compounds and 3830 transitions from SWATH 28 (batch 0 out of 1)
Thread 0_0 will analyze 766 compounds and 3830 transitions from SWATH 28 (batch 1 out of 1)
    87.88 %               Thread 0_0 will analyze 661 compounds and 3305 transitions from SWATH 29 (batch 0 out of 1)
Thread 0_0 will analyze 661 compounds and 3305 transitions from SWATH 29 (batch 1 out of 1)
    90.91 %               Thread 0_0 will analyze 603 compounds and 3015 transitions from SWATH 30 (batch 0 out of 1)
Thread 0_0 will analyze 603 compounds and 3015 transitions from SWATH 30 (batch 1 out of 1)
    93.94 %               Thread 0_0 will analyze 461 compounds and 2305 transitions from SWATH 31 (batch 0 out of 1)
Thread 0_0 will analyze 461 compounds and 2305 transitions from SWATH 31 (batch 1 out of 1)
    96.97 %               Thread 0_0 will analyze 351 compounds and 1755 transitions from SWATH 32 (batch 0 out of 1)
Thread 0_0 will analyze 351 compounds and 1755 transitions from SWATH 32 (batch 1 out of 1)
    100.00 %               
  -- done [took 16:38 m (CPU), 16:38 m (Wall)] -- 
OpenSwathWorkflow took 17:59 m (wall), 17:58 m (CPU), 33.37 s (system), 17:25 m (user); Peak Memory Usage: 680 MB.
Validate provided Swath windows file:
Read Swath window header: 'start	end'
Read Swath window file with 32 SWATH windows.
Progress of 'Load TSV file':
-- done [took 5.24 s (CPU), 5.26 s (Wall)] -- 
Loaded 8029 proteins, 95466 compounds with 477330 transitions.
Loading mzML file olgas_K121026_007_SW_Wayne_R2_d00-Mtb_Wayne_07.mzML using readoptions cache
Progress of 'Loading metadata file olgas_K121026_007_SW_Wayne_R2_d00-Mtb_Wayne_07.mzML':
Will analyze the metadata first to determine the number of SWATH windows and the window sizes.
Determined there to be 32 SWATH windows and in total 2277 MS1 spectra
Determined there to be 32 SWATH windows and in total 2277 MS1 spectra
-- done [took 14.35 s (CPU), 14.36 s (Wall)] -- 
Progress of 'Loading data file olgas_K121026_007_SW_Wayne_R2_d00-Mtb_Wayne_07.mzML':
Read chromatogram while reading SWATH files, did not expect that!

  -- done [took 01:00 m (CPU), 01:00 m (Wall)] -- 
Read Swath window header: 'start	end'
Read Swath window file with 32 SWATH windows.
Re-annotate from file: SWATH 400 / 425 (raw data) is annotated via swath_windows_file with 400 / 424.5
Re-annotate from file: SWATH 424 / 450 (raw data) is annotated via swath_windows_file with 424.5 / 449.5
Re-annotate from file: SWATH 449 / 475 (raw data) is annotated via swath_windows_file with 449.5 / 474.5
Re-annotate from file: SWATH 474 / 500 (raw data) is annotated via swath_windows_file with 474.5 / 499.5
Re-annotate from file: SWATH 499 / 525 (raw data) is annotated via swath_windows_file with 499.5 / 524.5
Re-annotate from file: SWATH 524 / 550 (raw data) is annotated via swath_windows_file with 524.5 / 549.5
Re-annotate from file: SWATH 549 / 575 (raw data) is annotated via swath_windows_file with 549.5 / 574.5
Re-annotate from file: SWATH 574 / 600 (raw data) is annotated via swath_windows_file with 574.5 / 599.5
Re-annotate from file: SWATH 599 / 625 (raw data) is annotated via swath_windows_file with 599.5 / 624.5
Re-annotate from file: SWATH 624 / 650 (raw data) is annotated via swath_windows_file with 624.5 / 649.5
Re-annotate from file: SWATH 649 / 675 (raw data) is annotated via swath_windows_file with 649.5 / 674.5
Re-annotate from file: SWATH 674 / 700 (raw data) is annotated via swath_windows_file with 674.5 / 699.5
Re-annotate from file: SWATH 699 / 725 (raw data) is annotated via swath_windows_file with 699.5 / 724.5
Re-annotate from file: SWATH 724 / 750 (raw data) is annotated via swath_windows_file with 724.5 / 749.5
Re-annotate from file: SWATH 749 / 775 (raw data) is annotated via swath_windows_file with 749.5 / 774.5
Re-annotate from file: SWATH 774 / 800 (raw data) is annotated via swath_windows_file with 774.5 / 799.5
Re-annotate from file: SWATH 799 / 825 (raw data) is annotated via swath_windows_file with 799.5 / 824.5
Re-annotate from file: SWATH 824 / 850 (raw data) is annotated via swath_windows_file with 824.5 / 849.5
Re-annotate from file: SWATH 849 / 875 (raw data) is annotated via swath_windows_file with 849.5 / 874.5
Re-annotate from file: SWATH 874 / 900 (raw data) is annotated via swath_windows_file with 874.5 / 899.5
Re-annotate from file: SWATH 899 / 925 (raw data) is annotated via swath_windows_file with 899.5 / 924.5
Re-annotate from file: SWATH 924 / 950 (raw data) is annotated via swath_windows_file with 924.5 / 949.5
Re-annotate from file: SWATH 949 / 975 (raw data) is annotated via swath_windows_file with 949.5 / 974.5
Re-annotate from file: SWATH 974 / 1000 (raw data) is annotated via swath_windows_file with 974.5 / 999.5
Re-annotate from file: SWATH 999 / 1025 (raw data) is annotated via swath_windows_file with 999.5 / 1024.5
Re-annotate from file: SWATH 1024 / 1050 (raw data) is annotated via swath_windows_file with 1024.5 / 1049.5
Re-annotate from file: SWATH 1049 / 1075 (raw data) is annotated via swath_windows_file with 1049.5 / 1074.5
Re-annotate from file: SWATH 1074 / 1100 (raw data) is annotated via swath_windows_file with 1074.5 / 1099.5
Re-annotate from file: SWATH 1099 / 1125 (raw data) is annotated via swath_windows_file with 1099.5 / 1124.5
Re-annotate from file: SWATH 1124 / 1150 (raw data) is annotated via swath_windows_file with 1124.5 / 1149.5
Re-annotate from file: SWATH 1149 / 1175 (raw data) is annotated via swath_windows_file with 1149.5 / 1174.5
Re-annotate from file: SWATH 1174 / 1200 (raw data) is annotated via swath_windows_file with 1174.5 / 1200
Will load iRT transitions and try to find iRT peptides

  Progress of 'Load TraML file':

  -- done [took 0.05 s (CPU), 0.05 s (Wall)] -- 

  Progress of 'Extract iRT chromatograms':

  -- done [took 0.37 s (CPU), 0.37 s (Wall)] -- 

  Progress of 'Retention time normalization':
Will analyse 20 peptides with a total of 100 transitions 
WARNING in SignalToNoiseEstimatorMedian: 1.84453% of all Signal-to-Noise estimates are too high, because the median was found in the rightmost histogram-bin. You should consider increasing 'max_intensity' (and maybe 'bin_count' with it, to keep bin width reasonable)
rsq: 0.999617 points: 20

  -- done [took 0.66 s (CPU), 0.66 s (Wall)] -- 
Will analyze 477330 transitions in total.

  Progress of 'Extracting and scoring transitions':
Use non-nested loop with 1 threads.
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 0 out of 3)
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 1 out of 3)
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 2 out of 3)
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 3 out of 3)
    6.06 %               Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 0 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 1 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 2 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 3 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 4 out of 4)
    9.09 %               Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 0 out of 5)
WARNING in SignalToNoiseEstimatorMedian: 100% of all windows were sparse. You should consider increasing 'win_len' or decreasing 'min_required_elements'
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 1 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 2 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 3 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 4 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 5 out of 5)
    12.12 %               Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 0 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 1 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 2 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 3 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 4 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 5 out of 5)
    15.15 %               Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 0 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 1 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 2 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 3 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 4 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 5 out of 5)
    18.18 %               Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 0 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 1 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 2 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 3 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 4 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 5 out of 5)
    21.21 %               Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 0 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 1 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 2 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 3 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 4 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 5 out of 5)
    24.24 %               Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 0 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 1 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 2 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 3 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 4 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 5 out of 5)
    27.27 %               Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 0 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 1 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 2 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 3 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 4 out of 4)
    30.30 %               Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 0 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 1 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 2 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 3 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 4 out of 4)
    33.33 %               Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 0 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 1 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 2 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 3 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 4 out of 4)
    36.36 %               Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 0 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 1 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 2 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 3 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 4 out of 4)
    39.39 %               Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 0 out of 3)
Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 1 out of 3)
Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 2 out of 3)
Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 3 out of 3)
    42.42 %               Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 0 out of 3)
Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 1 out of 3)
Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 2 out of 3)
Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 3 out of 3)
    45.45 %               Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 0 out of 3)
Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 1 out of 3)
Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 2 out of 3)
Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 3 out of 3)
    48.48 %               Thread 0_0 will analyze 2972 compounds and 14860 transitions from SWATH 16 (batch 0 out of 2)
Thread 0_0 will analyze 2972 compounds and 14860 transitions from SWATH 16 (batch 1 out of 2)
Thread 0_0 will analyze 2972 compounds and 14860 transitions from SWATH 16 (batch 2 out of 2)
    51.52 %               Thread 0_0 will analyze 2689 compounds and 13445 transitions from SWATH 17 (batch 0 out of 2)
Thread 0_0 will analyze 2689 compounds and 13445 transitions from SWATH 17 (batch 1 out of 2)
Thread 0_0 will analyze 2689 compounds and 13445 transitions from SWATH 17 (batch 2 out of 2)
    54.55 %               Thread 0_0 will analyze 2391 compounds and 11955 transitions from SWATH 18 (batch 0 out of 2)
Thread 0_0 will analyze 2391 compounds and 11955 transitions from SWATH 18 (batch 1 out of 2)
Thread 0_0 will analyze 2391 compounds and 11955 transitions from SWATH 18 (batch 2 out of 2)
    57.58 %               Thread 0_0 will analyze 2171 compounds and 10855 transitions from SWATH 19 (batch 0 out of 2)
Thread 0_0 will analyze 2171 compounds and 10855 transitions from SWATH 19 (batch 1 out of 2)
Thread 0_0 will analyze 2171 compounds and 10855 transitions from SWATH 19 (batch 2 out of 2)
    60.61 %               Thread 0_0 will analyze 2019 compounds and 10095 transitions from SWATH 20 (batch 0 out of 2)
Thread 0_0 will analyze 2019 compounds and 10095 transitions from SWATH 20 (batch 1 out of 2)
Thread 0_0 will analyze 2019 compounds and 10095 transitions from SWATH 20 (batch 2 out of 2)
    63.64 %               Thread 0_0 will analyze 1829 compounds and 9145 transitions from SWATH 21 (batch 0 out of 1)
Thread 0_0 will analyze 1829 compounds and 9145 transitions from SWATH 21 (batch 1 out of 1)
    66.67 %               Thread 0_0 will analyze 1628 compounds and 8140 transitions from SWATH 22 (batch 0 out of 1)
Thread 0_0 will analyze 1628 compounds and 8140 transitions from SWATH 22 (batch 1 out of 1)
    69.70 %               Thread 0_0 will analyze 1526 compounds and 7630 transitions from SWATH 23 (batch 0 out of 1)
Thread 0_0 will analyze 1526 compounds and 7630 transitions from SWATH 23 (batch 1 out of 1)
    72.73 %               Thread 0_0 will analyze 1371 compounds and 6855 transitions from SWATH 24 (batch 0 out of 1)
Thread 0_0 will analyze 1371 compounds and 6855 transitions from SWATH 24 (batch 1 out of 1)
    75.76 %               Thread 0_0 will analyze 1185 compounds and 5925 transitions from SWATH 25 (batch 0 out of 1)
Thread 0_0 will analyze 1185 compounds and 5925 transitions from SWATH 25 (batch 1 out of 1)
    78.79 %               Thread 0_0 will analyze 1106 compounds and 5530 transitions from SWATH 26 (batch 0 out of 1)
Thread 0_0 will analyze 1106 compounds and 5530 transitions from SWATH 26 (batch 1 out of 1)
    81.82 %               Thread 0_0 will analyze 889 compounds and 4445 transitions from SWATH 27 (batch 0 out of 1)
Thread 0_0 will analyze 889 compounds and 4445 transitions from SWATH 27 (batch 1 out of 1)
    84.85 %               Thread 0_0 will analyze 766 compounds and 3830 transitions from SWATH 28 (batch 0 out of 1)
Thread 0_0 will analyze 766 compounds and 3830 transitions from SWATH 28 (batch 1 out of 1)
    87.88 %               Thread 0_0 will analyze 661 compounds and 3305 transitions from SWATH 29 (batch 0 out of 1)
Thread 0_0 will analyze 661 compounds and 3305 transitions from SWATH 29 (batch 1 out of 1)
    90.91 %               Thread 0_0 will analyze 603 compounds and 3015 transitions from SWATH 30 (batch 0 out of 1)
Thread 0_0 will analyze 603 compounds and 3015 transitions from SWATH 30 (batch 1 out of 1)
    93.94 %               Thread 0_0 will analyze 461 compounds and 2305 transitions from SWATH 31 (batch 0 out of 1)
Thread 0_0 will analyze 461 compounds and 2305 transitions from SWATH 31 (batch 1 out of 1)
    96.97 %               Thread 0_0 will analyze 351 compounds and 1755 transitions from SWATH 32 (batch 0 out of 1)
Thread 0_0 will analyze 351 compounds and 1755 transitions from SWATH 32 (batch 1 out of 1)
    100.00 %               
  -- done [took 16:26 m (CPU), 16:27 m (Wall)] -- 
OpenSwathWorkflow took 17:48 m (wall), 17:48 m (CPU), 32.97 s (system), 17:15 m (user); Peak Memory Usage: 681 MB.
<WARNING in SignalToNoiseEstimatorMedian: 100% of all windows were sparse. You should consider increasing 'win_len' or decreasing 'min_required_elements'> occurred 10 times
Validate provided Swath windows file:
Read Swath window header: 'start	end'
Read Swath window file with 32 SWATH windows.
Progress of 'Load TSV file':
-- done [took 5.15 s (CPU), 5.14 s (Wall)] -- 
Loaded 8029 proteins, 95466 compounds with 477330 transitions.
Loading mzML file olgas_K121026_013_SW_Wayne_R3_d00-Mtb_Wayne_13.mzML using readoptions cache
Progress of 'Loading metadata file olgas_K121026_013_SW_Wayne_R3_d00-Mtb_Wayne_13.mzML':
Will analyze the metadata first to determine the number of SWATH windows and the window sizes.
Determined there to be 32 SWATH windows and in total 2277 MS1 spectra
Determined there to be 32 SWATH windows and in total 2277 MS1 spectra
-- done [took 14.23 s (CPU), 14.23 s (Wall)] -- 
Progress of 'Loading data file olgas_K121026_013_SW_Wayne_R3_d00-Mtb_Wayne_13.mzML':
Read chromatogram while reading SWATH files, did not expect that!

  -- done [took 59.65 s (CPU), 59.75 s (Wall)] -- 
Read Swath window header: 'start	end'
Read Swath window file with 32 SWATH windows.
Re-annotate from file: SWATH 400 / 425 (raw data) is annotated via swath_windows_file with 400 / 424.5
Re-annotate from file: SWATH 424 / 450 (raw data) is annotated via swath_windows_file with 424.5 / 449.5
Re-annotate from file: SWATH 449 / 475 (raw data) is annotated via swath_windows_file with 449.5 / 474.5
Re-annotate from file: SWATH 474 / 500 (raw data) is annotated via swath_windows_file with 474.5 / 499.5
Re-annotate from file: SWATH 499 / 525 (raw data) is annotated via swath_windows_file with 499.5 / 524.5
Re-annotate from file: SWATH 524 / 550 (raw data) is annotated via swath_windows_file with 524.5 / 549.5
Re-annotate from file: SWATH 549 / 575 (raw data) is annotated via swath_windows_file with 549.5 / 574.5
Re-annotate from file: SWATH 574 / 600 (raw data) is annotated via swath_windows_file with 574.5 / 599.5
Re-annotate from file: SWATH 599 / 625 (raw data) is annotated via swath_windows_file with 599.5 / 624.5
Re-annotate from file: SWATH 624 / 650 (raw data) is annotated via swath_windows_file with 624.5 / 649.5
Re-annotate from file: SWATH 649 / 675 (raw data) is annotated via swath_windows_file with 649.5 / 674.5
Re-annotate from file: SWATH 674 / 700 (raw data) is annotated via swath_windows_file with 674.5 / 699.5
Re-annotate from file: SWATH 699 / 725 (raw data) is annotated via swath_windows_file with 699.5 / 724.5
Re-annotate from file: SWATH 724 / 750 (raw data) is annotated via swath_windows_file with 724.5 / 749.5
Re-annotate from file: SWATH 749 / 775 (raw data) is annotated via swath_windows_file with 749.5 / 774.5
Re-annotate from file: SWATH 774 / 800 (raw data) is annotated via swath_windows_file with 774.5 / 799.5
Re-annotate from file: SWATH 799 / 825 (raw data) is annotated via swath_windows_file with 799.5 / 824.5
Re-annotate from file: SWATH 824 / 850 (raw data) is annotated via swath_windows_file with 824.5 / 849.5
Re-annotate from file: SWATH 849 / 875 (raw data) is annotated via swath_windows_file with 849.5 / 874.5
Re-annotate from file: SWATH 874 / 900 (raw data) is annotated via swath_windows_file with 874.5 / 899.5
Re-annotate from file: SWATH 899 / 925 (raw data) is annotated via swath_windows_file with 899.5 / 924.5
Re-annotate from file: SWATH 924 / 950 (raw data) is annotated via swath_windows_file with 924.5 / 949.5
Re-annotate from file: SWATH 949 / 975 (raw data) is annotated via swath_windows_file with 949.5 / 974.5
Re-annotate from file: SWATH 974 / 1000 (raw data) is annotated via swath_windows_file with 974.5 / 999.5
Re-annotate from file: SWATH 999 / 1025 (raw data) is annotated via swath_windows_file with 999.5 / 1024.5
Re-annotate from file: SWATH 1024 / 1050 (raw data) is annotated via swath_windows_file with 1024.5 / 1049.5
Re-annotate from file: SWATH 1049 / 1075 (raw data) is annotated via swath_windows_file with 1049.5 / 1074.5
Re-annotate from file: SWATH 1074 / 1100 (raw data) is annotated via swath_windows_file with 1074.5 / 1099.5
Re-annotate from file: SWATH 1099 / 1125 (raw data) is annotated via swath_windows_file with 1099.5 / 1124.5
Re-annotate from file: SWATH 1124 / 1150 (raw data) is annotated via swath_windows_file with 1124.5 / 1149.5
Re-annotate from file: SWATH 1149 / 1175 (raw data) is annotated via swath_windows_file with 1149.5 / 1174.5
Re-annotate from file: SWATH 1174 / 1200 (raw data) is annotated via swath_windows_file with 1174.5 / 1200
Will load iRT transitions and try to find iRT peptides

  Progress of 'Load TraML file':

  -- done [took 0.04 s (CPU), 0.04 s (Wall)] -- 

  Progress of 'Extract iRT chromatograms':

  -- done [took 0.36 s (CPU), 0.36 s (Wall)] -- 

  Progress of 'Retention time normalization':
Will analyse 20 peptides with a total of 100 transitions 
WARNING in SignalToNoiseEstimatorMedian: 2.15195% of all Signal-to-Noise estimates are too high, because the median was found in the rightmost histogram-bin. You should consider increasing 'max_intensity' (and maybe 'bin_count' with it, to keep bin width reasonable)
WARNING in SignalToNoiseEstimatorMedian: 1.05402% of all Signal-to-Noise estimates are too high, because the median was found in the rightmost histogram-bin. You should consider increasing 'max_intensity' (and maybe 'bin_count' with it, to keep bin width reasonable)
rsq: 0.93044 points: 20
rsq: 0.999368 points: 19

  -- done [took 0.65 s (CPU), 0.64 s (Wall)] -- 
Will analyze 477330 transitions in total.

  Progress of 'Extracting and scoring transitions':
Use non-nested loop with 1 threads.
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 0 out of 3)
WARNING in SignalToNoiseEstimatorMedian: 100% of all windows were sparse. You should consider increasing 'win_len' or decreasing 'min_required_elements'
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 1 out of 3)
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 2 out of 3)
Thread 0_0 will analyze 3828 compounds and 19140 transitions from SWATH 1 (batch 3 out of 3)
    6.06 %               Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 0 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 1 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 2 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 3 out of 4)
Thread 0_0 will analyze 4802 compounds and 24010 transitions from SWATH 2 (batch 4 out of 4)
    9.09 %               Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 0 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 1 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 2 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 3 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 4 out of 5)
Thread 0_0 will analyze 5564 compounds and 27820 transitions from SWATH 3 (batch 5 out of 5)
    12.12 %               Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 0 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 1 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 2 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 3 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 4 out of 5)
Thread 0_0 will analyze 5736 compounds and 28680 transitions from SWATH 4 (batch 5 out of 5)
    15.15 %               Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 0 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 1 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 2 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 3 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 4 out of 5)
Thread 0_0 will analyze 5816 compounds and 29080 transitions from SWATH 5 (batch 5 out of 5)
    18.18 %               Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 0 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 1 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 2 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 3 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 4 out of 5)
Thread 0_0 will analyze 5647 compounds and 28235 transitions from SWATH 6 (batch 5 out of 5)
    21.21 %               Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 0 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 1 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 2 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 3 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 4 out of 5)
Thread 0_0 will analyze 5328 compounds and 26640 transitions from SWATH 7 (batch 5 out of 5)
    24.24 %               Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 0 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 1 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 2 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 3 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 4 out of 5)
Thread 0_0 will analyze 5001 compounds and 25005 transitions from SWATH 8 (batch 5 out of 5)
    27.27 %               Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 0 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 1 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 2 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 3 out of 4)
Thread 0_0 will analyze 4767 compounds and 23835 transitions from SWATH 9 (batch 4 out of 4)
    30.30 %               Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 0 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 1 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 2 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 3 out of 4)
Thread 0_0 will analyze 4639 compounds and 23195 transitions from SWATH 10 (batch 4 out of 4)
    33.33 %               Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 0 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 1 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 2 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 3 out of 4)
Thread 0_0 will analyze 4313 compounds and 21565 transitions from SWATH 11 (batch 4 out of 4)
    36.36 %               Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 0 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 1 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 2 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 3 out of 4)
Thread 0_0 will analyze 4115 compounds and 20575 transitions from SWATH 12 (batch 4 out of 4)
    39.39 %               Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 0 out of 3)
Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 1 out of 3)
Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 2 out of 3)
Thread 0_0 will analyze 3827 compounds and 19135 transitions from SWATH 13 (batch 3 out of 3)
    42.42 %               Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 0 out of 3)
Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 1 out of 3)
Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 2 out of 3)
Thread 0_0 will analyze 3502 compounds and 17510 transitions from SWATH 14 (batch 3 out of 3)
    45.45 %               Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 0 out of 3)
Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 1 out of 3)
Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 2 out of 3)
Thread 0_0 will analyze 3312 compounds and 16560 transitions from SWATH 15 (batch 3 out of 3)
    48.48 %               Thread 0_0 will analyze 2972 compounds and 14860 transitions from SWATH 16 (batch 0 out of 2)
Thread 0_0 will analyze 2972 compounds and 14860 transitions from SWATH 16 (batch 1 out of 2)
Thread 0_0 will analyze 2972 compounds and 14860 transitions from SWATH 16 (batch 2 out of 2)
    51.52 %               Thread 0_0 will analyze 2689 compounds and 13445 transitions from SWATH 17 (batch 0 out of 2)
Thread 0_0 will analyze 2689 compounds and 13445 transitions from SWATH 17 (batch 1 out of 2)
Thread 0_0 will analyze 2689 compounds and 13445 transitions from SWATH 17 (batch 2 out of 2)
    54.55 %               Thread 0_0 will analyze 2391 compounds and 11955 transitions from SWATH 18 (batch 0 out of 2)
Thread 0_0 will analyze 2391 compounds and 11955 transitions from SWATH 18 (batch 1 out of 2)
Thread 0_0 will analyze 2391 compounds and 11955 transitions from SWATH 18 (batch 2 out of 2)
    57.58 %               Thread 0_0 will analyze 2171 compounds and 10855 transitions from SWATH 19 (batch 0 out of 2)
Thread 0_0 will analyze 2171 compounds and 10855 transitions from SWATH 19 (batch 1 out of 2)
Thread 0_0 will analyze 2171 compounds and 10855 transitions from SWATH 19 (batch 2 out of 2)
    60.61 %               Thread 0_0 will analyze 2019 compounds and 10095 transitions from SWATH 20 (batch 0 out of 2)
Thread 0_0 will analyze 2019 compounds and 10095 transitions from SWATH 20 (batch 1 out of 2)
Thread 0_0 will analyze 2019 compounds and 10095 transitions from SWATH 20 (batch 2 out of 2)
    63.64 %               Thread 0_0 will analyze 1829 compounds and 9145 transitions from SWATH 21 (batch 0 out of 1)
Thread 0_0 will analyze 1829 compounds and 9145 transitions from SWATH 21 (batch 1 out of 1)
    66.67 %               Thread 0_0 will analyze 1628 compounds and 8140 transitions from SWATH 22 (batch 0 out of 1)
Thread 0_0 will analyze 1628 compounds and 8140 transitions from SWATH 22 (batch 1 out of 1)
    69.70 %               Thread 0_0 will analyze 1526 compounds and 7630 transitions from SWATH 23 (batch 0 out of 1)
Thread 0_0 will analyze 1526 compounds and 7630 transitions from SWATH 23 (batch 1 out of 1)
    72.73 %               Thread 0_0 will analyze 1371 compounds and 6855 transitions from SWATH 24 (batch 0 out of 1)
Thread 0_0 will analyze 1371 compounds and 6855 transitions from SWATH 24 (batch 1 out of 1)
    75.76 %               Thread 0_0 will analyze 1185 compounds and 5925 transitions from SWATH 25 (batch 0 out of 1)
Thread 0_0 will analyze 1185 compounds and 5925 transitions from SWATH 25 (batch 1 out of 1)
    78.79 %               Thread 0_0 will analyze 1106 compounds and 5530 transitions from SWATH 26 (batch 0 out of 1)
Thread 0_0 will analyze 1106 compounds and 5530 transitions from SWATH 26 (batch 1 out of 1)
    81.82 %               Thread 0_0 will analyze 889 compounds and 4445 transitions from SWATH 27 (batch 0 out of 1)
Thread 0_0 will analyze 889 compounds and 4445 transitions from SWATH 27 (batch 1 out of 1)
    84.85 %               Thread 0_0 will analyze 766 compounds and 3830 transitions from SWATH 28 (batch 0 out of 1)
Thread 0_0 will analyze 766 compounds and 3830 transitions from SWATH 28 (batch 1 out of 1)
    87.88 %               Thread 0_0 will analyze 661 compounds and 3305 transitions from SWATH 29 (batch 0 out of 1)
Thread 0_0 will analyze 661 compounds and 3305 transitions from SWATH 29 (batch 1 out of 1)
    90.91 %               Thread 0_0 will analyze 603 compounds and 3015 transitions from SWATH 30 (batch 0 out of 1)
Thread 0_0 will analyze 603 compounds and 3015 transitions from SWATH 30 (batch 1 out of 1)
    93.94 %               Thread 0_0 will analyze 461 compounds and 2305 transitions from SWATH 31 (batch 0 out of 1)
Thread 0_0 will analyze 461 compounds and 2305 transitions from SWATH 31 (batch 1 out of 1)
    96.97 %               Thread 0_0 will analyze 351 compounds and 1755 transitions from SWATH 32 (batch 0 out of 1)
Thread 0_0 will analyze 351 compounds and 1755 transitions from SWATH 32 (batch 1 out of 1)
    100.00 %               
  -- done [took 16:50 m (CPU), 16:55 m (Wall)] -- 
OpenSwathWorkflow took 18:16 m (wall), 18:11 m (CPU), 35.68 s (system), 17:35 m (user); Peak Memory Usage: 684 MB.
<WARNING in SignalToNoiseEstimatorMedian: 100% of all windows were sparse. You should consider increasing 'win_len' or decreasing 'min_required_elements'> occurred 25 times
