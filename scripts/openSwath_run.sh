#!/bin/bash

#STR = "ile"
#echo STR
echo "Run OpenSwath."
for filename in *.mzML
		do
				STR=${filename:16}
				if [ ${STR::5} == "32fix" ]; then 
						echo $filename
						OpenSwathWorkflow -in $filename -tr ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_fixed_curated_target_decoy.pqp -out_osw 20210207_osw_output.$filename.osw  -tempDirectory tmp -readOptions cacheWorkingInMemory -batchSize 1000 -Scoring:stop_report_after_feature -1 -min_upper_edge_dist 1 -extra_rt_extraction_window 100 -min_rsq 0.95 -min_coverage 0.6 -Scoring:Scores:use_dia_scores true -rt_extraction_window 600 -mz_extraction_window 30 -threads 6
				fi
				if [ ${STR::5} == "32var" ]; then
						echo $filename
						OpenSwathWorkflow -in $filename -tr ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_var_curated_target_decoy.pqp -out_osw 20210207_osw_output.$filename.osw  -tempDirectory tmp -readOptions cacheWorkingInMemory -batchSize 1000 -Scoring:stop_report_after_feature -1 -min_upper_edge_dist 1 -extra_rt_extraction_window 100 -min_rsq 0.95 -min_coverage 0.6 -Scoring:Scores:use_dia_scores true -rt_extraction_window 600 -mz_extraction_window 30 -threads 6
				fi
				if [ ${STR::5} == "64fix" ]; then
						echo $filename
						OpenSwathWorkflow -in $filename -tr ecolihumanyeast_concat_mayu_IRR_cons_openswath_64w_fixed_curated_target_decoy.pqp -out_osw 20210207_osw_output.$filename.osw  -tempDirectory tmp -readOptions cacheWorkingInMemory -batchSize 1000 -Scoring:stop_report_after_feature -1 -min_upper_edge_dist 1 -extra_rt_extraction_window 100 -min_rsq 0.95 -min_coverage 0.6 -Scoring:Scores:use_dia_scores true -rt_extraction_window 600 -mz_extraction_window 30 -threads 6
				fi
				if [ ${STR::5} == "64var" ]; then
						echo $filename
						OpenSwathWorkflow -in $filename -tr ecolihumanyeast_concat_mayu_IRR_cons_openswath_64w_var_curated_target_decoy.pqp -out_osw 20210207_osw_output.$filename.osw  -tempDirectory tmp -readOptions cacheWorkingInMemory -batchSize 1000 -Scoring:stop_report_after_feature -1 -min_upper_edge_dist 1 -extra_rt_extraction_window 100 -min_rsq 0.95 -min_coverage 0.6 -Scoring:Scores:use_dia_scores true -rt_extraction_window 600 -mz_extraction_window 30 -threads 6
				fi

				#echo ${STR::5}
				#OpenSwathWorkflow -in $filename -tr ecolihumanyeast_concat_mayu_IRR_cons_openswath_32w_var_curated_target_decoy.pqp -out_osw target_decoy_osw_output.$filename.osw  -tempDirectory tmp -readOptions cacheWorkingInMemory -batchSize 1000 -Scoring:stop_report_after_feature -1 -min_upper_edge_dist 1 -extra_rt_extraction_window 100 -min_rsq 0.95 -min_coverage 0.6 -Scoring:Scores:use_dia_scores true -rt_extraction_window 600 -mz_extraction_window 30 -threads 6
				#echo "Generated ${filename::-6}_target_decoy.TraML"
		done;


