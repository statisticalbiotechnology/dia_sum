#!/bin/bash

echo "Statistical validation using Percolator..."
for filename in 20210207*
		do			
				PercolatorAdapter -in_osw $filename -out percolator_output_$filename.osw 
				echo "Generated percolator_output_${filename}"
		done;


