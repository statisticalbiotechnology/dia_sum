#!/bin/bash

echo "Converting Percolator processed .osw file to .tsv..."
for filename in percolator_output*
		do			
				pyprophet export --in=${filename}				 
				echo "Generated ${filename::-4}.tsv"
		done;


