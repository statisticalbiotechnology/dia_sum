
(py27) ptruong@planck:/hdd_14T/data/PASS00779/ftp.peptideatlas.org/experiment_20200211$ ./pyprophet_run.sh 
Run pyprophet

used parameters:

    apply_scorer                               : None
    apply_weights                              : None
    compute.probabilities                      : False
    d_score.cutoff                             : 0.5
    delim.in                                   : '\t'
    delim.out                                  : '\t'
    export.mayu                                : False
    final_statistics.fdr_all_pg                : False
    final_statistics.lambda                    : 0.4
    ignore.invalid_score_columns               : 1
    is_test                                    : 0
    num_processes                              : 1
    semi_supervised_learner.initial_fdr        : 0.15
    semi_supervised_learner.initial_lambda     : 0.4
    semi_supervised_learner.iteration_fdr      : 0.02
    semi_supervised_learner.iteration_lambda   : 0.4
    semi_supervised_learner.num_iter           : 5
    target.dir                                 : None
    target.overwrite                           : 0
    xeval.fraction                             : 0.5
    xeval.num_iter                             : 5

INFO -- [pid=401041] : 2021-02-11 23:14:19,152: config settings:
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     apply_scorer: None
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     apply_weights: None
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     compute.probabilities: False
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     d_score.cutoff: 0.5
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     delim.in: 	
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     delim.out: 	
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     export.mayu: False
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     final_statistics.fdr_all_pg: False
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     final_statistics.lambda: 0.4
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     ignore.invalid_score_columns: 1
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     is_test: 0
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     num_processes: 1
INFO -- [pid=401041] : 2021-02-11 23:14:19,152:     semi_supervised_learner.initial_fdr: 0.15
INFO -- [pid=401041] : 2021-02-11 23:14:19,153:     semi_supervised_learner.initial_lambda: 0.4
INFO -- [pid=401041] : 2021-02-11 23:14:19,153:     semi_supervised_learner.iteration_fdr: 0.02
INFO -- [pid=401041] : 2021-02-11 23:14:19,153:     semi_supervised_learner.iteration_lambda: 0.4
INFO -- [pid=401041] : 2021-02-11 23:14:19,153:     semi_supervised_learner.num_iter: 5
INFO -- [pid=401041] : 2021-02-11 23:14:19,153:     target.dir: None
INFO -- [pid=401041] : 2021-02-11 23:14:19,153:     target.overwrite: 0
INFO -- [pid=401041] : 2021-02-11 23:14:19,153:     xeval.fraction: 0.5
INFO -- [pid=401041] : 2021-02-11 23:14:19,153:     xeval.num_iter: 5
INFO -- [pid=401041] : 2021-02-11 23:14:19,153: read osw_001.tsv
INFO -- [pid=401041] : 2021-02-11 23:14:28,088: learn and apply scorer to osw_001.tsv
WARNING -- [pid=401041] : 2021-02-11 23:14:28,231: %s. pyprophet skips this.
WARNING -- [pid=401041] : 2021-02-11 23:14:28,234: %s. pyprophet skips this.
WARNING -- [pid=401041] : 2021-02-11 23:14:28,237: %s. pyprophet skips this.
WARNING -- [pid=401041] : 2021-02-11 23:14:28,240: %s. pyprophet skips this.
WARNING -- [pid=401041] : 2021-02-11 23:14:28,247: %s. pyprophet skips this.
INFO -- [pid=401041] : 2021-02-11 23:14:29,329: data set contains 46979 decoy and 47830 target transition groups
INFO -- [pid=401041] : 2021-02-11 23:14:29,343: summary input file:
INFO -- [pid=401041] : 2021-02-11 23:14:29,343:    473847 lines
INFO -- [pid=401041] : 2021-02-11 23:14:29,372:    94809 transition groups
INFO -- [pid=401041] : 2021-02-11 23:14:29,372:    27 scores including main score
INFO -- [pid=401041] : 2021-02-11 23:14:29,373: start 5 cross evals using 1 processes
INFO -- [pid=401041] : 2021-02-11 23:14:29,373: start learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:31,119: end learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:31,122: start learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:32,031: end learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:32,033: start learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:32,904: end learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:32,908: start learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:33,848: end learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:33,850: start learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:34,732: end learn_randomized
INFO -- [pid=401041] : 2021-02-11 23:14:34,736: finished cross evals
INFO -- [pid=401041] : 2021-02-11 23:14:34,968: mean m_score = 2.921214e-01, std_dev m_score = 1.109264e-01
INFO -- [pid=401041] : 2021-02-11 23:14:34,969: mean s_value = 9.584961e-01, std_dev s_value = 1.459572e-01
INFO -- [pid=401041] : 2021-02-11 23:14:35,433: calculated scoring and statistics
INFO -- [pid=401041] : 2021-02-11 23:14:35,438: processing osw_001.tsv finished
INFO -- [pid=401041] : 2021-02-11 23:14:35,438: time needed: 00:00:16.3

==============================================================================

   qvalue    svalue       TP       FP    ...              FN       FDR      sens    cutoff
0    0.00  0.000040      2.0      0.0    ...     29713.90923  0.000000  0.000067  8.157814
1    0.01  0.729613  21681.0    219.0    ...      8034.90923  0.010000  0.729609  2.251974
2    0.02  0.773922  22998.0    469.0    ...      6717.90923  0.019993  0.773929  1.943555
3    0.05  0.848305  25208.0   1327.0    ...      4507.90923  0.050008  0.848300  1.437669
4    0.10  0.918814  27303.0   3034.0    ...      2412.90923  0.099997  0.918801  0.952013
5    0.20  1.000000  29802.0   7450.0    ...       -86.09077  0.199999  1.000000  0.219745
6    0.30  1.000000  30688.0  13152.0    ...      -972.09077  0.299991  1.000000 -0.598941
7    0.40  1.000000  29716.0  18114.0    ...        -0.09077  0.378718  1.000000 -5.341124
8    0.50       NaN      NaN      NaN    ...             NaN       NaN       NaN       NaN

[9 rows x 9 columns]

==============================================================================

WRITTEN:  osw_001_summary_stat.csv
WRITTEN:  osw_001_full_stat.csv
WRITTEN:  osw_001_report.pdf
WRITTEN:  osw_001_cutoffs.txt
WRITTEN:  osw_001_svalues.txt
WRITTEN:  osw_001_qvalues.txt
WRITTEN:  osw_001_dscores_top_target_peaks.txt
WRITTEN:  osw_001_dscores_top_decoy_peaks.txt
WRITTEN:  osw_001_with_dscore.csv
WRITTEN:  osw_001_with_dscore_filtered.csv
WRITTEN:  osw_001_scorer.bin
WRITTEN:  osw_001_weights.txt

NEEDED 16 seconds and 321 msecs wall time


used parameters:

    apply_scorer                               : None
    apply_weights                              : None
    compute.probabilities                      : False
    d_score.cutoff                             : 0.5
    delim.in                                   : '\t'
    delim.out                                  : '\t'
    export.mayu                                : False
    final_statistics.fdr_all_pg                : False
    final_statistics.lambda                    : 0.4
    ignore.invalid_score_columns               : 1
    is_test                                    : 0
    num_processes                              : 1
    semi_supervised_learner.initial_fdr        : 0.15
    semi_supervised_learner.initial_lambda     : 0.4
    semi_supervised_learner.iteration_fdr      : 0.02
    semi_supervised_learner.iteration_lambda   : 0.4
    semi_supervised_learner.num_iter           : 5
    target.dir                                 : None
    target.overwrite                           : 0
    xeval.fraction                             : 0.5
    xeval.num_iter                             : 5

INFO -- [pid=401061] : 2021-02-11 23:15:03,072: config settings:
INFO -- [pid=401061] : 2021-02-11 23:15:03,072:     apply_scorer: None
INFO -- [pid=401061] : 2021-02-11 23:15:03,072:     apply_weights: None
INFO -- [pid=401061] : 2021-02-11 23:15:03,072:     compute.probabilities: False
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     d_score.cutoff: 0.5
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     delim.in: 	
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     delim.out: 	
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     export.mayu: False
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     final_statistics.fdr_all_pg: False
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     final_statistics.lambda: 0.4
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     ignore.invalid_score_columns: 1
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     is_test: 0
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     num_processes: 1
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     semi_supervised_learner.initial_fdr: 0.15
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     semi_supervised_learner.initial_lambda: 0.4
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     semi_supervised_learner.iteration_fdr: 0.02
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     semi_supervised_learner.iteration_lambda: 0.4
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     semi_supervised_learner.num_iter: 5
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     target.dir: None
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     target.overwrite: 0
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     xeval.fraction: 0.5
INFO -- [pid=401061] : 2021-02-11 23:15:03,073:     xeval.num_iter: 5
INFO -- [pid=401061] : 2021-02-11 23:15:03,073: read osw_007.tsv
INFO -- [pid=401061] : 2021-02-11 23:15:15,080: learn and apply scorer to osw_007.tsv
WARNING -- [pid=401061] : 2021-02-11 23:15:15,226: %s. pyprophet skips this.
WARNING -- [pid=401061] : 2021-02-11 23:15:15,230: %s. pyprophet skips this.
WARNING -- [pid=401061] : 2021-02-11 23:15:15,233: %s. pyprophet skips this.
WARNING -- [pid=401061] : 2021-02-11 23:15:15,236: %s. pyprophet skips this.
WARNING -- [pid=401061] : 2021-02-11 23:15:15,243: %s. pyprophet skips this.
INFO -- [pid=401061] : 2021-02-11 23:15:16,233: data set contains 46976 decoy and 47828 target transition groups
INFO -- [pid=401061] : 2021-02-11 23:15:16,248: summary input file:
INFO -- [pid=401061] : 2021-02-11 23:15:16,249:    473680 lines
INFO -- [pid=401061] : 2021-02-11 23:15:16,276:    94804 transition groups
INFO -- [pid=401061] : 2021-02-11 23:15:16,276:    27 scores including main score
INFO -- [pid=401061] : 2021-02-11 23:15:16,276: start 5 cross evals using 1 processes
INFO -- [pid=401061] : 2021-02-11 23:15:16,276: start learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:17,909: end learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:17,912: start learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:18,862: end learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:18,865: start learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:19,820: end learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:19,823: start learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:20,786: end learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:20,790: start learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:21,764: end learn_randomized
INFO -- [pid=401061] : 2021-02-11 23:15:21,768: finished cross evals
INFO -- [pid=401061] : 2021-02-11 23:15:21,998: mean m_score = 2.740820e-01, std_dev m_score = 1.060573e-01
INFO -- [pid=401061] : 2021-02-11 23:15:21,999: mean s_value = 9.579899e-01, std_dev s_value = 1.475656e-01
INFO -- [pid=401061] : 2021-02-11 23:15:22,470: calculated scoring and statistics
INFO -- [pid=401061] : 2021-02-11 23:15:22,475: processing osw_007.tsv finished
INFO -- [pid=401061] : 2021-02-11 23:15:22,475: time needed: 00:00:19.4

==============================================================================

   qvalue    svalue       TP       FP    ...              FN           FDR      sens    cutoff
0    0.00  0.000039      1.0      0.0    ...     30745.61905  6.321116e-12  0.000033  8.095052
1    0.01  0.755414  23226.0    235.0    ...      7520.61905  9.998016e-03  0.755400  2.199747
2    0.02  0.798072  24538.0    501.0    ...      6208.61905  2.001261e-02  0.798071  1.890091
3    0.05  0.865202  26602.0   1400.0    ...      4144.61905  4.999433e-02  0.865201  1.392231
4    0.10  0.932505  28671.0   3186.0    ...      2075.61905  9.999563e-02  0.932493  0.894397
5    0.20  1.000000  31103.0   7775.0    ...      -356.38095  1.999866e-01  1.000000  0.122831
6    0.30  1.000000  31648.0  13564.0    ...      -901.38095  2.999995e-01  1.000000 -0.816741
7    0.40  1.000000  30747.0  17081.0    ...        -0.38095  3.571419e-01  1.000000 -9.018684
8    0.50       NaN      NaN      NaN    ...             NaN           NaN       NaN       NaN

[9 rows x 9 columns]

==============================================================================

WRITTEN:  osw_007_summary_stat.csv
WRITTEN:  osw_007_full_stat.csv
WRITTEN:  osw_007_report.pdf
WRITTEN:  osw_007_cutoffs.txt
WRITTEN:  osw_007_svalues.txt
WRITTEN:  osw_007_qvalues.txt
WRITTEN:  osw_007_dscores_top_target_peaks.txt
WRITTEN:  osw_007_dscores_top_decoy_peaks.txt
WRITTEN:  osw_007_with_dscore.csv
WRITTEN:  osw_007_with_dscore_filtered.csv
WRITTEN:  osw_007_scorer.bin
WRITTEN:  osw_007_weights.txt

NEEDED 19 seconds and 436 msecs wall time


used parameters:

    apply_scorer                               : None
    apply_weights                              : None
    compute.probabilities                      : False
    d_score.cutoff                             : 0.5
    delim.in                                   : '\t'
    delim.out                                  : '\t'
    export.mayu                                : False
    final_statistics.fdr_all_pg                : False
    final_statistics.lambda                    : 0.4
    ignore.invalid_score_columns               : 1
    is_test                                    : 0
    num_processes                              : 1
    semi_supervised_learner.initial_fdr        : 0.15
    semi_supervised_learner.initial_lambda     : 0.4
    semi_supervised_learner.iteration_fdr      : 0.02
    semi_supervised_learner.iteration_lambda   : 0.4
    semi_supervised_learner.num_iter           : 5
    target.dir                                 : None
    target.overwrite                           : 0
    xeval.fraction                             : 0.5
    xeval.num_iter                             : 5

INFO -- [pid=401075] : 2021-02-11 23:15:48,034: config settings:
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     apply_scorer: None
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     apply_weights: None
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     compute.probabilities: False
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     d_score.cutoff: 0.5
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     delim.in: 	
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     delim.out: 	
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     export.mayu: False
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     final_statistics.fdr_all_pg: False
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     final_statistics.lambda: 0.4
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     ignore.invalid_score_columns: 1
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     is_test: 0
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     num_processes: 1
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     semi_supervised_learner.initial_fdr: 0.15
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     semi_supervised_learner.initial_lambda: 0.4
INFO -- [pid=401075] : 2021-02-11 23:15:48,034:     semi_supervised_learner.iteration_fdr: 0.02
INFO -- [pid=401075] : 2021-02-11 23:15:48,035:     semi_supervised_learner.iteration_lambda: 0.4
INFO -- [pid=401075] : 2021-02-11 23:15:48,035:     semi_supervised_learner.num_iter: 5
INFO -- [pid=401075] : 2021-02-11 23:15:48,035:     target.dir: None
INFO -- [pid=401075] : 2021-02-11 23:15:48,035:     target.overwrite: 0
INFO -- [pid=401075] : 2021-02-11 23:15:48,035:     xeval.fraction: 0.5
INFO -- [pid=401075] : 2021-02-11 23:15:48,035:     xeval.num_iter: 5
INFO -- [pid=401075] : 2021-02-11 23:15:48,035: read osw_013.tsv
INFO -- [pid=401075] : 2021-02-11 23:15:56,632: learn and apply scorer to osw_013.tsv
WARNING -- [pid=401075] : 2021-02-11 23:15:56,777: %s. pyprophet skips this.
WARNING -- [pid=401075] : 2021-02-11 23:15:56,780: %s. pyprophet skips this.
WARNING -- [pid=401075] : 2021-02-11 23:15:56,783: %s. pyprophet skips this.
WARNING -- [pid=401075] : 2021-02-11 23:15:56,786: %s. pyprophet skips this.
WARNING -- [pid=401075] : 2021-02-11 23:15:56,793: %s. pyprophet skips this.
INFO -- [pid=401075] : 2021-02-11 23:15:57,793: data set contains 46969 decoy and 47819 target transition groups
INFO -- [pid=401075] : 2021-02-11 23:15:57,808: summary input file:
INFO -- [pid=401075] : 2021-02-11 23:15:57,808:    473604 lines
INFO -- [pid=401075] : 2021-02-11 23:15:57,833:    94788 transition groups
INFO -- [pid=401075] : 2021-02-11 23:15:57,833:    27 scores including main score
INFO -- [pid=401075] : 2021-02-11 23:15:57,833: start 5 cross evals using 1 processes
INFO -- [pid=401075] : 2021-02-11 23:15:57,833: start learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:15:58,968: end learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:15:58,971: start learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:15:59,929: end learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:15:59,931: start learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:16:00,889: end learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:16:00,891: start learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:16:01,845: end learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:16:01,847: start learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:16:02,816: end learn_randomized
INFO -- [pid=401075] : 2021-02-11 23:16:02,818: finished cross evals
INFO -- [pid=401075] : 2021-02-11 23:16:03,046: mean m_score = 3.089100e-01, std_dev m_score = 1.153751e-01
INFO -- [pid=401075] : 2021-02-11 23:16:03,047: mean s_value = 9.596269e-01, std_dev s_value = 1.438004e-01
INFO -- [pid=401075] : 2021-02-11 23:16:03,517: calculated scoring and statistics
INFO -- [pid=401075] : 2021-02-11 23:16:03,521: processing osw_013.tsv finished
INFO -- [pid=401075] : 2021-02-11 23:16:03,522: time needed: 00:00:15.5

==============================================================================

   qvalue    svalue       TP       FP    ...               FN           FDR      sens    cutoff
0    0.00  0.000070      2.0      0.0    ...     28731.934223  6.356471e-12  0.000070  8.004934
1    0.01  0.725878  20857.0    211.0    ...      7876.934223  1.000101e-02  0.725866  2.299484
2    0.02  0.771711  22174.0    453.0    ...      6559.934223  2.000220e-02  0.771701  1.992645
3    0.05  0.840140  24141.0   1270.0    ...      4592.934223  4.999490e-02  0.840156  1.512056
4    0.10  0.908083  26093.0   2899.0    ...      2640.934223  9.999381e-02  0.908090  1.029004
5    0.20  0.994929  28588.0   7147.0    ...       145.934223  2.000024e-01  0.994921  0.323526
6    0.30  1.000000  29841.0  12789.0    ...     -1107.065777  3.000058e-01  1.000000 -0.438010
7    0.40  1.000000  28734.0  19085.0    ...        -0.065777  3.991105e-01  1.000000 -7.562698
8    0.50       NaN      NaN      NaN    ...              NaN           NaN       NaN       NaN

[9 rows x 9 columns]

==============================================================================

WRITTEN:  osw_013_summary_stat.csv
WRITTEN:  osw_013_full_stat.csv
WRITTEN:  osw_013_report.pdf
WRITTEN:  osw_013_cutoffs.txt
WRITTEN:  osw_013_svalues.txt
WRITTEN:  osw_013_qvalues.txt
WRITTEN:  osw_013_dscores_top_target_peaks.txt
WRITTEN:  osw_013_dscores_top_decoy_peaks.txt
WRITTEN:  osw_013_with_dscore.csv
WRITTEN:  osw_013_with_dscore_filtered.csv
WRITTEN:  osw_013_scorer.bin
WRITTEN:  osw_013_weights.txt

NEEDED 15 seconds and 519 msecs wall time
